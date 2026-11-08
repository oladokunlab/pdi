#!/usr/bin/env Rscript
# =============================================================================
# 19_sparse_pls.R
#
# Sparse PLS: finds co-varying modules between the whole genus matrix and the
# whole gut-health domain matrix at once, rather than one genus and one
# collapsed score at a time.
#
# RUN:  Rscript 19_sparse_pls.R     (needs 04b and 06 to have run)
# TIME: 5-30 minutes (the cross-validation step is slow).
#
# -----------------------------------------------------------------------------
# WHY THIS IS A DIFFERENT QUESTION FROM EVERYTHING ELSE IN THIS PROJECT
# -----------------------------------------------------------------------------
# MaAsLin2 and LASSO both regress the microbiome against ONE outcome column
# (GHI_Combined). Dysbiosis is a community property, not necessarily reducible
# to a single number - a set of genera might co-vary strongly with a
# PARTICULAR combination of domains (say, barrier and inflammation together)
# without any single genus surviving a test against the combined score. sPLS
# looks for that directly: latent components that explain shared variation
# between the ENTIRE genus block and the ENTIRE domain block simultaneously.
#
# -----------------------------------------------------------------------------
# THE Y BLOCK, AND WHY IT IS NOT SIMPLY THE 10 RAW DOMAIN SCORES
# -----------------------------------------------------------------------------
# The source template uses five raw, always-available metrics as Y. This
# project's GHI has ten domains, most of which are scored on only a subset of
# samples (SCFA n=74-75, Immunoglobulin_Systemic n=63-84), because the
# underlying biomarker panel is sparse in exactly the way that forced the
# domain/proration structure in the first place - complete cases across every
# domain is zero, the same problem the registered PCA method could not survive.
#
# sPLS itself requires a complete-case Y block (no missingness). Rather than
# assume ten domains jointly have enough complete cases, this script requires
# it and PROGRESSIVELY DROPS the domain with the fewest complete cases,
# refitting the sample count each time, until either a workable Y block is
# found or too few domains remain to be informative. Every drop is logged, so
# which domains ended up in the final Y block - and why - is visible rather
# than assumed.
#
# -----------------------------------------------------------------------------
# UNIT OF ANALYSIS
# -----------------------------------------------------------------------------
# Arm-level CLR abundances (the same X used throughout this pipeline), matched
# to arm-level domain scores. Consistent with every other stage; the source
# template operates on individual samples, which would repeat the
# pseudo-replication problem already fixed elsewhere in this project.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("mixOmics", quietly = TRUE))
  stop("Package 'mixOmics' required.\n",
       "  if (!require(\"BiocManager\")) install.packages(\"BiocManager\")\n",
       "  BiocManager::install(\"mixOmics\", lib = Sys.getenv(\"R_LIBS_USER\"))")
suppressPackageStartupMessages(library(mixOmics))

MIN_COMPLETE <- 30    # floor below which a Y block is not attempted
NCOMP <- 3
KEEPX_GRID <- c(5, 10, 15, 20, 30, 50)
set.seed(42)

log_msg("[1/6] Loading X (genus CLR) and Y (domain scores)")
mf <- file.path(R_OUT_DIR, "gdi", "gdi_model.rds")
if (!file.exists(mf)) stop("Not found: ", mf, "\n  Run 06_lasso_selection.R first.")
m <- readRDS(mf)
X_all <- m$X; arm_ids <- rownames(X_all)
log_msg("  X: %d arms x %d genera", nrow(X_all), ncol(X_all))

gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
dom_cols <- grep("^GHI_dom_", names(gh), value = TRUE)
if (!length(dom_cols)) stop("No GHI_dom_* columns found in ghi_scores_pca.tsv")
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
Y_all <- sapply(dom_cols, function(cn) {
  v <- suppressWarnings(as.numeric(gh[[cn]]))
  unname(tapply(v, gk, mean, na.rm = TRUE)[arm_ids])
})
rownames(Y_all) <- arm_ids
colnames(Y_all) <- sub("^GHI_dom_", "", dom_cols)
log_msg("  Y: %d domains available (%s)", ncol(Y_all), paste(colnames(Y_all), collapse = ", "))

# =============================================================================
# 2 - BUILD A COMPLETE-CASE Y BLOCK, DROPPING THE SPARSEST DOMAIN AS NEEDED
# =============================================================================
log_msg("[2/6] Finding a usable complete-case Y block")
active <- colnames(Y_all)
repeat {
  ok <- stats::complete.cases(Y_all[, active, drop = FALSE])
  n_ok <- sum(ok)
  log_msg("  domains: %s -> %d complete-case arms", paste(active, collapse = ", "), n_ok)
  if (n_ok >= MIN_COMPLETE || length(active) <= 2) break
  miss <- colSums(is.na(Y_all[, active, drop = FALSE]))
  drop_dom <- names(which.max(miss))
  log_msg("  dropping '%s' (most missing: %d of %d) and retrying", drop_dom, miss[drop_dom], nrow(Y_all))
  active <- setdiff(active, drop_dom)
}
if (n_ok < MIN_COMPLETE) {
  log_msg("  *** Could not reach %d complete cases even with only %d domains (%s).",
          MIN_COMPLETE, length(active), paste(active, collapse = ", "))
  log_msg("  *** This is itself informative: the biomarker panel is too sparse to")
  log_msg("  *** support a multivariate Y block at all, independent of sPLS's")
  log_msg("  *** ability to find structure in it.")
  quit(status = 0)
}
common <- arm_ids[ok]
X <- X_all[common, , drop = FALSE]
Y <- Y_all[common, active, drop = FALSE]
log_msg("  final: %d arms, %d genera (X), %d domains (Y): %s",
        nrow(X), ncol(X), ncol(Y), paste(active, collapse = ", "))

# Drop zero-variance genus columns, which tune.spls cannot handle
keep_x <- apply(X, 2, function(v) stats::sd(v) > 0)
X <- X[, keep_x, drop = FALSE]
log_msg("  %d genera retained after removing zero-variance columns", ncol(X))

# =============================================================================
# 3 - TUNE SPARSITY VIA CROSS-VALIDATION
# =============================================================================
log_msg("[3/6] Tuning sparsity (5-fold CV, 10 repeats) - this is the slow step")
ncomp_use <- min(NCOMP, ncol(Y) - 1, floor(nrow(X) / 10))
if (ncomp_use < 1) ncomp_use <- 1
log_msg("  using %d component(s) given %d arms and %d Y columns", ncomp_use, nrow(X), ncol(Y))

keepx_grid <- KEEPX_GRID[KEEPX_GRID <= ncol(X)]
tune_res <- tryCatch(
  mixOmics::tune.spls(X, Y, ncomp = ncomp_use, test.keepX = keepx_grid,
                      validation = "Mfold", folds = min(5, nrow(X) %/% 6),
                      nrepeat = 10, progressBar = FALSE),
  error = function(e) { log_msg("  tune.spls failed: %s", conditionMessage(e)); NULL })

chosen_keepx <- if (!is.null(tune_res)) tune_res$choice.keepX else rep(min(keepx_grid), ncomp_use)
log_msg("  chosen keepX per component: %s", paste(chosen_keepx, collapse = ", "))

# =============================================================================
# 4 - FINAL MODEL
# =============================================================================
log_msg("[4/6] Fitting final sPLS model")
keepy_use <- rep(min(ncol(Y), 5), ncomp_use)
spls_fit <- mixOmics::spls(X, Y, ncomp = ncomp_use, keepX = chosen_keepx, keepY = keepy_use)

# =============================================================================
# 5 - REPORT LOADINGS AND CROSS-CHECK AGAINST OTHER METHODS
# =============================================================================
log_msg("[5/6] Top genera per component")
gdi_stab <- if (file.exists(file.path(R_OUT_DIR, "gdi", "gdi_stability.tsv")))
  fread(file.path(R_OUT_DIR, "gdi", "gdi_stability.tsv")) else NULL

for (cc in seq_len(ncomp_use)) {
  lx <- spls_fit$loadings$X[, cc]
  lx <- lx[lx != 0]
  lx <- sort(lx, decreasing = TRUE)
  ly <- spls_fit$loadings$Y[, cc]; ly <- ly[ly != 0]
  log_msg("")
  log_msg("  --- Component %d ---", cc)
  log_msg("  domain loadings: %s",
          paste(sprintf("%s=%+.2f", names(ly), ly), collapse = "  "))
  log_msg("  top genera (by |loading|):")
  top <- head(lx[order(-abs(lx))], 10)
  for (g in names(top)) {
    stab <- if (!is.null(gdi_stab)) {
      r <- gdi_stab[genus == g]
      if (nrow(r)) sprintf("  [LASSO stability %.0f%%]", 100 * r$stability[1]) else ""
    } else ""
    log_msg("    %-32s loading=%+.3f%s", gsub("^g__", "", g), top[g], stab)
  }
}

# =============================================================================
# 6 - WRITE
# =============================================================================
log_msg("[6/6] Writing")
out_dir <- file.path(R_OUT_DIR, "gdi", "spls")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(spls_fit, file.path(out_dir, "spls_fit.rds"))

lines <- c("SPARSE PLS: MICROBIOME vs GUT HEALTH DOMAINS", strrep("=", 74), "",
  sprintf("Arms: %d, genera: %d, domains in Y: %s",
          nrow(X), ncol(X), paste(active, collapse = ", ")),
  sprintf("Components: %d, keepX: %s", ncomp_use, paste(chosen_keepx, collapse = ", ")), "")
for (cc in seq_len(ncomp_use)) {
  lx <- spls_fit$loadings$X[, cc]; lx <- lx[lx != 0]
  ly <- spls_fit$loadings$Y[, cc]; ly <- ly[ly != 0]
  lines <- c(lines, sprintf("--- Component %d ---", cc),
            sprintf("domains: %s", paste(sprintf("%s=%+.2f", names(ly), ly), collapse = "  ")),
            sprintf("genera: %s",
                    paste(sprintf("%s=%+.3f", gsub("^g__", "", names(lx)), lx), collapse = "  ")), "")
}
writeLines(lines, file.path(out_dir, "spls_summary.txt"))
log_msg("summary -> %s", file.path(out_dir, "spls_summary.txt"))
log_msg("fit object -> %s", file.path(out_dir, "spls_fit.rds"))
log_msg("=== DONE ===")
