#!/usr/bin/env Rscript
# =============================================================================
# 06_lasso_selection.R
#
# Stage 2 of GDI construction: from the genera MaAsLin2 flagged, select a
# sparse, stable panel by penalised regression with bootstrap stability
# selection.
#
# RUN:  Rscript 06_lasso_selection.R      (needs 05_maaslin2.R first)
# TIME: a few minutes.
#
# -----------------------------------------------------------------------------
# THE UNIT OF ANALYSIS: WHY THIS FITS AT ARM LEVEL, NOT RUN LEVEL
# -----------------------------------------------------------------------------
# The outcome (GHI) is a published group mean. Every sequencing run within a
# treatment arm carries the SAME GHI value - the papers report one number per
# arm, not one per bird. Fitting at run level therefore adds no information
# about the outcome while multiplying the apparent sample size.
#
# In this corpus that distortion is severe. Biagi_2020 contributes 714 runs
# (23% of all runs) but only a handful of distinct GHI values; Miska_2024
# contributes 360. A run-level fit would let those two studies dominate the
# panel purely through replication, and every p-value and cross-validation
# estimate would be computed against an n that does not exist.
#
# Abundances are therefore averaged to one row per study x treatment arm - the
# resolution at which the outcome actually varies (171 arms across the corpus).
# This is the standard remedy for pseudo-replication and it makes the
# effective sample size honest, at the cost of appearing much smaller. The
# smaller number is the real one.
#
# -----------------------------------------------------------------------------
# WHY BOOTSTRAP RESAMPLING IS BY STUDY
# -----------------------------------------------------------------------------
# Arms within a study share husbandry, diet base, bird genetics, sequencing
# run and laboratory protocol, so they are not independent. Resampling arms
# would treat them as though they were, understating uncertainty and
# overstating how stable the selected panel is.
#
# Stability selection therefore subsamples STUDIES (50% without replacement,
# per Meinshausen and Buhlmann), refits, and records which genera survive.
# A genus counts as stable if it is selected in at least the threshold share
# of subsamples - meaning it holds up when half the evidence base is removed.
#
# -----------------------------------------------------------------------------
# COMPOSITIONALITY
# -----------------------------------------------------------------------------
# Inputs are CLR-transformed (03_batch_correction.R). CLR is the log-ratio
# transform that makes ordinary regression valid on compositional data: a
# coefficient set summing to zero on CLR values is algebraically identical to
# Aitchison's linear log-contrast model. glmnet does not impose that zero-sum
# constraint, so the fit is an approximation to the log-contrast model rather
# than exactly it. The zeroSum or FLORAL packages impose it directly and are
# the more rigorous option if available; the approximation is noted here
# rather than left implicit.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

# ---- registered parameters ---------------------------------------------------
N_BOOT          <- 1000   # registered
STABILITY_THR   <- 0.70   # registered; 60-90% reported as sensitivity
SUBSAMPLE_FRAC  <- 0.50   # Meinshausen-Buhlmann standard
LAMBDA_RULE     <- "1se"  # registered primary; "min" reported alongside
set.seed(42)              # reproducibility

log_msg("[1/6] Loading MaAsLin2 output")
inp <- file.path(R_OUT_DIR, "gdi", "gdi_input.rds")
if (!file.exists(inp)) stop(sprintf("Not found: %s\nRun 05_maaslin2.R first.", inp))
d <- readRDS(inp)
clr <- d$clr; meta <- d$meta; selected <- d$selected
log_msg("  %d genera x %d runs; %d flagged by MaAsLin2", nrow(clr), ncol(clr), length(selected))
if (length(selected) < 2) stop("Fewer than 2 genera passed MaAsLin2 - nothing to select from.")

# ---- restrict to the MaAsLin2 survivors --------------------------------------
keep <- rownames(clr) %in% selected
clr <- clr[keep, , drop = FALSE]
log_msg("  restricted to %d genera", nrow(clr))

# =============================================================================
# STEP 2 - COLLAPSE TO ONE ROW PER STUDY x TREATMENT ARM
# =============================================================================
log_msg("[2/6] Collapsing runs to study x arm (the resolution of the outcome)")
if (!"Treatment_Group" %in% names(meta)) {
  # 05 stores Study_ID/Segment/GHI; recover the arm from the crosswalk
  cw <- as.data.frame(read_excel(METADATA_XLSX, sheet = "Crosswalk_16S"),
                      stringsAsFactors = FALSE)
  meta$Treatment_Group <- cw$Treatment_Group[match(rownames(meta), trimws(cw$Run))]
}
arm <- paste(meta$Study_ID, meta$Treatment_Group, sep = "||")
ok <- !is.na(arm) & !is.na(meta$GHI)
clr <- clr[, ok, drop = FALSE]; meta <- meta[ok, , drop = FALSE]; arm <- arm[ok]

ua <- unique(arm)
X <- t(sapply(ua, function(a) rowMeans(clr[, arm == a, drop = FALSE], na.rm = TRUE)))
y <- vapply(ua, function(a) mean(meta$GHI[arm == a]), numeric(1))
study <- vapply(ua, function(a) meta$Study_ID[arm == a][1], character(1))
n_runs <- vapply(ua, function(a) sum(arm == a), numeric(1))

log_msg("  %d runs -> %d study-arm rows across %d studies",
        ncol(clr), nrow(X), length(unique(study)))
log_msg("  runs per arm: median %.0f (range %d-%d)",
        stats::median(n_runs), min(n_runs), max(n_runs))
log_msg("  effective sample size for the model is %d, not %d", nrow(X), ncol(clr))

drop <- !stats::complete.cases(X) | is.na(y)
if (any(drop)) {
  log_msg("  dropping %d arm(s) with incomplete data", sum(drop))
  X <- X[!drop, , drop = FALSE]; y <- y[!drop]; study <- study[!drop]
}
if (nrow(X) < 20) log_msg("  *** WARNING: only %d arms - selection will be unstable", nrow(X))

# =============================================================================
# STEP 3 - LASSO ON THE FULL DATA
# =============================================================================
log_msg("[3/6] LASSO (cross-validated)")

# -----------------------------------------------------------------------------
# FOLD CONSTRUCTION - explicit, not left to cv.glmnet's default
#
# cv.glmnet's default fold assignment is a simple random split. With few rows
# (n=158 here) and a continuous but not-fully-distinct outcome, a random split
# can by chance put many identical or near-identical y values in one fold,
# collapsing that fold's variance towards zero and destabilising the CV error
# curve - the previous run of this script produced lambda.min == lambda.1se
# and a CV/null error ratio of 1.009, both symptoms of exactly this.
#
# Folds are constructed here by sorting on y and assigning stratified blocks,
# which spreads the outcome range evenly across folds and removes this failure
# mode. A minimum fold count is also enforced so no fold is too small to
# estimate a variance from at all.
# -----------------------------------------------------------------------------
n <- nrow(X)
n_folds <- max(3, min(10, floor(n / 8)))   # at least ~8 rows per fold
ord <- order(y)
fold_id <- integer(n)
fold_id[ord] <- ((seq_len(n) - 1) %% n_folds) + 1
log_msg("  using %d folds (stratified on y, ~%.0f rows/fold) instead of a random split",
        n_folds, n / n_folds)

fold_var <- tapply(y, fold_id, stats::var)
if (any(fold_var == 0 | is.na(fold_var))) {
  log_msg("  *** WARNING: %d fold(s) still have zero variance in y even after",
          sum(fold_var == 0 | is.na(fold_var)))
  log_msg("  *** stratification. Reducing fold count further.")
  n_folds <- max(3, n_folds - 2)
  fold_id[ord] <- ((seq_len(n) - 1) %% n_folds) + 1
}

cvf <- cv.glmnet(X, y, alpha = 1, foldid = fold_id, standardize = TRUE)
lam <- if (LAMBDA_RULE == "1se") cvf$lambda.1se else cvf$lambda.min

if (isTRUE(all.equal(cvf$lambda.min, cvf$lambda.1se))) {
  log_msg("  *** lambda.min equals lambda.1se (%.5f) - the CV curve is flat.", cvf$lambda.min)
  log_msg("  *** This means cross-validated error does not distinguish any")
  log_msg("  *** amount of shrinkage from the null model at this fold structure.")
}

get_coefs <- function(fit, l) {
  co <- as.matrix(coef(fit, s = l))
  co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  co[co[, 1] != 0, 1]
}
sel_1se <- get_coefs(cvf, cvf$lambda.1se)
sel_min <- get_coefs(cvf, cvf$lambda.min)
log_msg("  lambda.1se selects %d genera (registered primary)", length(sel_1se))
log_msg("  lambda.min selects %d genera (sensitivity)", length(sel_min))

# =============================================================================
# STEP 4 - BOOTSTRAP STABILITY SELECTION, SUBSAMPLING STUDIES
# =============================================================================
log_msg("[4/6] Stability selection: %d iterations, %.0f%% of studies each",
        N_BOOT, 100 * SUBSAMPLE_FRAC)
studies <- unique(study)
n_take <- max(3, floor(length(studies) * SUBSAMPLE_FRAC))
counts <- setNames(numeric(ncol(X)), colnames(X))
n_ok <- 0
for (b in seq_len(N_BOOT)) {
  take <- sample(studies, n_take)
  idx <- which(study %in% take)
  if (length(idx) < 10 || length(unique(y[idx])) < 5) next
  fit <- tryCatch(
    cv.glmnet(X[idx, , drop = FALSE], y[idx], alpha = 1,
              nfolds = min(5, length(idx)), standardize = TRUE),
    error = function(e) NULL)
  if (is.null(fit)) next
  s <- get_coefs(fit, if (LAMBDA_RULE == "1se") fit$lambda.1se else fit$lambda.min)
  counts[names(s)] <- counts[names(s)] + 1
  n_ok <- n_ok + 1
}
freq <- counts / max(n_ok, 1)
log_msg("  %d of %d iterations usable", n_ok, N_BOOT)

stable <- sort(freq[freq >= STABILITY_THR], decreasing = TRUE)
log_msg("  %d genera stable at >= %.0f%%", length(stable), 100 * STABILITY_THR)

# Threshold sensitivity, reported so the 70% cut is visibly not arbitrary
log_msg("")
log_msg("  panel size across thresholds:")
for (th in c(0.6, 0.7, 0.8, 0.9)) {
  log_msg("    >= %.0f%%  %2d genera", 100 * th, sum(freq >= th))
}

log_msg("")
log_msg("  stable panel (selection frequency, LASSO coefficient):")
if (length(stable) == 0) {
  log_msg("    NONE reach %.0f%%. The strongest are:", 100 * STABILITY_THR)
  for (g in names(sort(freq, decreasing = TRUE))[1:min(10, length(freq))])
    log_msg("      %-42s %.0f%%", g, 100 * freq[[g]])
} else {
  full <- get_coefs(cvf, lam)
  for (g in names(stable)) {
    cf <- if (g %in% names(full)) full[[g]] else NA_real_
    log_msg("    %-42s %3.0f%%  coef=%s", g, 100 * stable[[g]],
            if (is.na(cf)) "  (not at this lambda)" else sprintf("%+.4f", cf))
  }
}

# =============================================================================
# STEP 5 - IN-SAMPLE FIT
#
# Reported only as a floor check. If the model cannot predict its own training
# data the target is noise; a good in-sample fit proves nothing about
# generalisation, which is what 08_loso_validate.R measures.
# =============================================================================
log_msg("[5/6] In-sample fit (a floor check, not evidence of generalisation)")
pred <- as.numeric(predict(cvf, newx = X, s = lam))
r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
log_msg("  R2=%.3f  Pearson r=%.3f  Spearman rho=%.3f",
        r2, stats::cor(y, pred), stats::cor(y, pred, method = "spearman"))
if (r2 < 0.05) {
  log_msg("  *** R2 below 0.05 in-sample. The genera carry almost no information")
  log_msg("  *** about this target even before generalisation is tested.")
}

# =============================================================================
# STEP 6 - EXPORT
# =============================================================================
log_msg("[6/6] Writing outputs")
panel <- if (length(stable) > 0) names(stable) else names(sort(freq, decreasing = TRUE))[1:min(10, length(freq))]
gd <- file.path(R_OUT_DIR, "gdi"); dir.create(gd, recursive = TRUE, showWarnings = FALSE)
writeLines(panel, file.path(gd, "gdi_panel.txt"))
fwrite(data.table(genus = names(freq), stability = as.numeric(freq),
                  coef_1se = as.numeric(sapply(names(freq), function(g)
                    if (g %in% names(sel_1se)) sel_1se[[g]] else 0))
       )[order(-stability)], file.path(gd, "gdi_stability.tsv"), sep = "\t")
saveRDS(list(X = X, y = y, study = study, panel = panel, freq = freq,
             lambda = lam, cvfit = cvf, target = d$target),
        file.path(gd, "gdi_model.rds"))

out <- c("=== LASSO taxon selection with bootstrap stability ===", "",
  sprintf("Target: %s", d$target),
  sprintf("Genera entering from MaAsLin2: %d", length(selected)),
  sprintf("Study-arm rows: %d across %d studies", nrow(X), length(unique(study))),
  sprintf("Runs collapsed: %d -> %d rows", ncol(clr), nrow(X)),
  "",
  sprintf("lambda.1se selects %d; lambda.min selects %d", length(sel_1se), length(sel_min)),
  sprintf("Stable at >=%.0f%%: %d genera", 100*STABILITY_THR, length(stable)),
  "",
  "Panel size by threshold:",
  paste0("  >=", c(60,70,80,90), "%: ", sapply(c(.6,.7,.8,.9), function(t) sum(freq>=t)), " genera"),
  "",
  sprintf("In-sample R2 = %.3f (floor check only)", r2),
  "",
  "--- method notes ---",
  "Fitted at study x arm resolution because the outcome is a published group",
  "mean: every run in an arm shares one GHI value, so run-level fitting adds",
  "no information while inflating n. Biagi_2020 alone contributes 714 runs;",
  "without collapsing, it would dominate selection through replication.",
  "",
  "Stability selection subsamples STUDIES, not arms, because arms within a",
  "study share husbandry, genetics, sequencing run and protocol. A genus is",
  "stable if it survives when half the evidence base is removed.",
  "",
  "CLR-transformed input makes ordinary regression valid on compositional",
  "data. glmnet does not impose the zero-sum coefficient constraint, so this",
  "approximates rather than exactly implements Aitchison's log-contrast model;",
  "zeroSum or FLORAL would impose it directly.",
  "",
  "In-sample R2 is a floor check. Generalisation is measured only by",
  "leave-one-study-out validation in 08_loso_validate.R.")
writeLines(out, file.path(gd, "lasso_summary.txt"))
log_msg("  panel     -> %s", file.path(gd, "gdi_panel.txt"))
log_msg("  stability -> %s", file.path(gd, "gdi_stability.tsv"))
log_msg("=== DONE ===")

# =============================================================================
# DIAGNOSTIC APPENDED FOR TROUBLESHOOTING A ZERO-GENE SELECTION
# Not part of the normal pipeline output. Prints the full lambda path so a
# null selection can be distinguished from a genuinely uninformative panel.
# =============================================================================
log_msg("")
log_msg("=== DIAGNOSTIC: lambda path detail ===")
log_msg("  lambda.min = %.5f  (%d genera at this lambda)", cvf$lambda.min,
        sum(get_coefs(cvf, cvf$lambda.min) != 0))
log_msg("  lambda.1se = %.5f  (%d genera at this lambda)", cvf$lambda.1se,
        sum(get_coefs(cvf, cvf$lambda.1se) != 0))
log_msg("  lambda.min CV error = %.5f", min(cvf$cvm))
log_msg("  null-model CV error (predict the mean) = %.5f", stats::var(y))
log_msg("  ratio (min CV error / null error): %.3f  (<1.0 means the model beats the mean)",
        min(cvf$cvm) / stats::var(y))
# how many genera survive at progressively looser lambda values
for (frac in c(1, 0.5, 0.2, 0.1, 0.05, 0.01)) {
  l <- cvf$lambda.min * frac
  n_sel <- sum(get_coefs(cvf, l) != 0)
  log_msg("  at lambda = %.5f (%.0f%% of lambda.min): %d genera selected", l, 100*frac, n_sel)
}

log_msg("")
log_msg("=== DIAGNOSTIC 2: is X degenerate? ===")
log_msg("  X dim: %d x %d", nrow(X), ncol(X))
log_msg("  any NA in X: %s", any(is.na(X)))
log_msg("  any NaN in X: %s", any(is.nan(X)))
log_msg("  any Inf in X: %s", any(is.infinite(X)))
col_sd <- apply(X, 2, stats::sd, na.rm = TRUE)
log_msg("  columns with zero/NA SD: %d of %d", sum(col_sd == 0 | is.na(col_sd)), ncol(X))
log_msg("  X value range: [%.4f, %.4f]", min(X, na.rm=TRUE), max(X, na.rm=TRUE))
log_msg("  y value range: [%.4f, %.4f]  y SD: %.4f", min(y), max(y), stats::sd(y))
log_msg("  first 3x3 corner of X:")
print(round(X[1:min(3,nrow(X)), 1:min(3,ncol(X))], 4))
log_msg("  first 5 y values: %s", paste(round(y[1:5], 3), collapse=", "))
