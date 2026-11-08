#!/usr/bin/env Rscript
# =============================================================================
# 19b_sparse_pls_fallback.R
#
# Sparse PLS using the CRAN 'spls' package instead of mixOmics.
#
# RUN:  Rscript 19b_sparse_pls_fallback.R    (needs 04b and 06 to have run)
#
# -----------------------------------------------------------------------------
# WHY A SECOND VERSION
# -----------------------------------------------------------------------------
# 19_sparse_pls.R (mixOmics) failed to install: RSpectra, one of mixOmics'
# dependencies, requires a C++ compiler flag (-std=gnu23) that this cluster's
# available GCC does not support. That is a toolchain limitation, not a code
# problem, and may not be fixable without an HPRC-side compiler upgrade.
#
# 'spls' (Chun & Keles, 2010) implements the same underlying method - sparse
# partial least squares - with a much lighter dependency chain (no Bioconductor,
# minimal compiled code), and is far more likely to install cleanly here.
#
# -----------------------------------------------------------------------------
# WHAT'S THE SAME AND WHAT'S DIFFERENT FROM 19_sparse_pls.R
# -----------------------------------------------------------------------------
# Same X (arm-level CLR genus matrix), same Y-block construction (domain
# scores, sparsest domains dropped until a workable complete-case block is
# found), same unit of analysis (arms, not runs).
#
# Different tuning mechanics: spls's cv.spls() tunes two parameters - eta
# (sparsity, 0 to 1: higher means fewer genera kept) and K (number of
# components) - by k-fold cross-validation, rather than mixOmics' keepX grid.
# The reported output (which genera load on which components, and how those
# components relate to the domains) answers the same question either way.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("spls", quietly = TRUE))
  stop("Package 'spls' required.\n  install.packages(\"spls\", lib = Sys.getenv(\"R_LIBS_USER\"), repos = \"https://cloud.r-project.org\")")
suppressPackageStartupMessages(library(spls))

MIN_COMPLETE <- 30
ETA_GRID <- seq(0.1, 0.9, by = 0.1)
K_GRID   <- 1:3
set.seed(42)

log_msg("[1/5] Loading X (genus CLR) and Y (domain scores)")
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

log_msg("[2/5] Finding a usable complete-case Y block")
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
  log_msg("  *** Could not reach %d complete cases. Biomarker panel too sparse", MIN_COMPLETE)
  log_msg("  *** for a multivariate Y block, independent of which sPLS package is used.")
  quit(status = 0)
}
common <- arm_ids[ok]
X <- X_all[common, , drop = FALSE]
Y <- Y_all[common, active, drop = FALSE]
keep_x <- apply(X, 2, function(v) stats::sd(v) > 0)
X <- X[, keep_x, drop = FALSE]
log_msg("  final: %d arms, %d genera, %d domains: %s",
        nrow(X), ncol(X), ncol(Y), paste(active, collapse = ", "))

log_msg("[3/5] Cross-validating eta (sparsity) and K (components)")
cv <- tryCatch(
  spls::cv.spls(X, Y, eta = ETA_GRID, K = K_GRID, plot.it = FALSE),
  error = function(e) { log_msg("  cv.spls failed: %s", conditionMessage(e)); NULL })
if (is.null(cv)) stop("Cross-validation failed; cannot proceed.")
log_msg("  chosen eta = %.1f, K = %d", cv$eta.opt, cv$K.opt)

log_msg("[4/5] Fitting final sPLS model")
fit <- spls::spls(X, Y, eta = cv$eta.opt, K = cv$K.opt)
print(spls::print.spls(fit))

sel_idx <- fit$A   # indices of selected (non-zero) predictors
sel_genera <- colnames(X)[sel_idx]
log_msg("  %d of %d genera retained at this sparsity level", length(sel_genera), ncol(X))

coefs <- coef(fit)  # genera x domains coefficient matrix
log_msg("")
log_msg("  Genera selected, by domain loading:")
for (g in sel_genera) {
  row <- coefs[g, ]
  nz <- row[row != 0]
  if (!length(nz)) next
  log_msg("    %-32s %s", gsub("^g__", "", g),
          paste(sprintf("%s=%+.3f", names(nz), nz), collapse = "  "))
}

log_msg("[5/5] Writing")
out_dir <- file.path(R_OUT_DIR, "gdi", "spls_fallback")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path(out_dir, "spls_fit.rds"))
out <- c("SPARSE PLS (spls package): MICROBIOME vs GUT HEALTH DOMAINS", strrep("=", 74), "",
  sprintf("Arms: %d, genera: %d, domains: %s", nrow(X), ncol(X), paste(active, collapse=", ")),
  sprintf("eta (sparsity) = %.1f, K (components) = %d", cv$eta.opt, cv$K.opt),
  sprintf("Genera selected: %d of %d", length(sel_genera), ncol(X)), "",
  "Genus coefficients by domain:")
for (g in sel_genera) {
  row <- coefs[g, ]; nz <- row[row != 0]
  if (length(nz)) out <- c(out, sprintf("  %-32s %s", gsub("^g__","",g),
                                       paste(sprintf("%s=%+.3f", names(nz), nz), collapse="  ")))
}
writeLines(out, file.path(out_dir, "spls_summary.txt"))
log_msg("summary -> %s", file.path(out_dir, "spls_summary.txt"))
log_msg("=== DONE ===")
