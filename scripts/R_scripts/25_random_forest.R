#!/usr/bin/env Rscript
# =============================================================================
# 25_random_forest.R
#
# Random forest instead of penalized linear regression, following the CRC
# meta-analysis approach (Pekel et al., Cell Host & Microbe 2026). Same arms,
# same target, same within-study centering, same study-level resampling.
#
# RUN:  Rscript 25_random_forest.R    (needs 06_lasso_selection.R; needs ranger)
# TIME: ~10 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS IS GENUINELY DIFFERENT FROM THE ALPHA SWEEP
# -----------------------------------------------------------------------------
# Elastic Net across alpha 1.00 -> 0.05 raised stability 31% -> 40% but left
# R2 at 0.000 throughout. That rules out one explanation (correlated genera
# competing under L1) but leaves another untested: every model tried so far,
# at any alpha, is LINEAR and ADDITIVE. It can only ever fit
# "GHI = b1*genus1 + b2*genus2 + ...".
#
# Random forest is neither. It can represent thresholds ("this genus only
# matters once it exceeds some abundance") and interactions ("high Blautia
# helps only when Enterococcus is low"). If the microbiome-GHI relationship
# has that shape, no linear model at any penalty would find it, and the flat
# R2 across the whole alpha sweep would look exactly as observed.
#
# The CRC meta-analysis used random forest for this reason and reported
# AUC >= 0.78 transferring between 16S and shotgun profiles - strong
# cross-study generalisation on a heterogeneous corpus.
#
# -----------------------------------------------------------------------------
# THE OVERFITTING RISK, AND WHY LOSO IS THE ACTUAL TEST HERE
# -----------------------------------------------------------------------------
# The CRC study had 6,779 samples. This corpus has 141 arms. Random forest is
# MORE prone to overfitting at small n than LASSO, not less - it will happily
# memorise 141 rows and report a high in-sample R2 that means nothing.
#
# So in-sample R2 is reported only as a diagnostic, and the number that
# decides the question is LEAVE-ONE-STUDY-OUT: train on all studies but one,
# predict the held-out study, pool predictions. If in-sample R2 is high and
# LOSO is near zero, that is overfitting and should be read as a negative
# result, not a positive one. That comparison is printed explicitly.
#
# Permutation importance (not impurity importance) is used, since impurity
# importance is biased toward high-cardinality continuous predictors.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("ranger", quietly = TRUE))
  stop("Package 'ranger' required.\n  install.packages(\"ranger\", lib = Sys.getenv(\"R_LIBS_USER\"), repos = \"https://cloud.r-project.org\")")
suppressPackageStartupMessages(library(ranger))

N_TREES <- 1000
N_BOOT  <- 500        # fewer than LASSO's 1000: each RF fit is far slower
SUBFRAC <- 0.50
STABILITY_THR <- 0.70
IMP_QUANTILE  <- 0.75  # a genus counts as "selected" in a resample if its
                       # importance is in the top quartile for that fit
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")

log_msg("[1/5] Loading the same arm-level data used by every other GDI script")
m <- readRDS(file.path(R_OUT_DIR, "gdi", "gdi_model.rds"))
X <- m$X; y <- m$y; study <- m$study
log_msg("  %d arms x %d genera across %d studies", nrow(X), ncol(X), length(unique(study)))

log_msg("[2/5] Within-study centering (identical to the registered method)")
Xc <- X; yc <- y
for (s in unique(study)) { i <- study == s
  Xc[i, ] <- sweep(X[i, , drop = FALSE], 2, colMeans(X[i, , drop = FALSE]), "-")
  yc[i] <- y[i] - mean(y[i]) }
df <- as.data.frame(Xc); names(df) <- make.names(colnames(Xc)); df$.y <- yc
name_map <- setNames(colnames(Xc), make.names(colnames(Xc)))

log_msg("[3/5] In-sample fit (diagnostic only - see header on overfitting)")
rf_full <- ranger(dependent.variable.name = ".y", data = df, num.trees = N_TREES,
                  importance = "permutation", seed = 42)
pred_in <- rf_full$predictions   # OOB predictions, already a mild guard
r2_oob <- max(0, 1 - sum((yc - pred_in)^2) / sum((yc - mean(yc))^2))
log_msg("  out-of-bag R2 = %.3f   (LASSO comparison: 0.040; Elastic Net: 0.000)", r2_oob)
log_msg("  OOB is already an honest-ish estimate, but LOSO below is the real test.")

log_msg("[4/5] Leave-one-study-out - the number that decides this")
studies <- unique(study)
pred_loso <- rep(NA_real_, length(yc))
for (s in studies) {
  tr <- study != s; te <- study == s
  if (sum(tr) < 20 || sum(te) < 1) next
  fit <- tryCatch(ranger(dependent.variable.name = ".y", data = df[tr, , drop = FALSE],
                         num.trees = N_TREES, seed = 42),
                 error = function(e) NULL)
  if (is.null(fit)) next
  pred_loso[te] <- predict(fit, data = df[te, , drop = FALSE])$predictions
}
ok <- !is.na(pred_loso)
r2_loso <- if (sum(ok) > 5) max(0, 1 - sum((yc[ok]-pred_loso[ok])^2)/sum((yc[ok]-mean(yc[ok]))^2)) else NA
rho_loso <- if (sum(ok) > 5) suppressWarnings(cor(yc[ok], pred_loso[ok], method="spearman")) else NA
log_msg("  LOSO R2  = %.3f  (on %d of %d arms)", r2_loso, sum(ok), length(yc))
log_msg("  LOSO rho = %+.3f", rho_loso)
if (!is.na(r2_loso) && !is.na(r2_oob)) {
  if (r2_oob > 0.15 && r2_loso < 0.05)
    log_msg("  *** OOB high, LOSO near zero: this is overfitting, NOT signal.")
  else if (r2_loso >= 0.05)
    log_msg("  *** LOSO R2 clears 0.05 - the first configuration in this project to do so.")
  else
    log_msg("  Both low: random forest finds no more than the linear models did.")
}

log_msg("[5/5] Stability selection on permutation importance, %d resamples", N_BOOT)
ntake <- max(3, floor(length(studies) * SUBFRAC))
cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
for (b in seq_len(N_BOOT)) {
  tk <- sample(studies, ntake); idx <- which(study %in% tk)
  if (length(idx) < 20) next
  fit <- tryCatch(ranger(dependent.variable.name = ".y", data = df[idx, , drop = FALSE],
                         num.trees = 300, importance = "permutation", seed = b),
                 error = function(e) NULL)
  if (is.null(fit)) next
  imp <- ranger::importance(fit)
  if (!length(imp) || all(!is.finite(imp))) next
  thr <- stats::quantile(imp, IMP_QUANTILE, na.rm = TRUE)
  sel <- name_map[names(imp)[imp >= thr & imp > 0]]
  sel <- sel[!is.na(sel)]
  cnt[sel] <- cnt[sel] + 1; nok <- nok + 1
}
freq <- if (nok > 0) cnt / nok else cnt
freq <- freq[!(names(freq) %in% UNASSIGNABLE)]
log_msg("  %d of %d resamples usable", nok, N_BOOT)
log_msg("  genera above the registered %.0f%% threshold: %d", 100*STABILITY_THR, sum(freq >= STABILITY_THR))
log_msg("  NOTE: this threshold is not directly comparable to the LASSO one - here")
log_msg("  it means 'in the top quartile by importance', not 'given a non-zero coefficient'.")
top <- sort(freq, decreasing = TRUE)[1:min(10, length(freq))]
for (i in seq_along(top)) log_msg("    %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i])

log_msg("")
log_msg("Overall permutation importance from the full model:")
imp_full <- sort(ranger::importance(rf_full), decreasing = TRUE)
for (i in seq_len(min(10, length(imp_full))))
  log_msg("    %-32s %.5f", gsub("^g__","",name_map[names(imp_full)[i]]), imp_full[i])

out <- c("RANDOM FOREST (non-linear) vs PENALIZED LINEAR MODELS", strrep("=", 74), "",
  sprintf("Arms: %d, genera: %d, studies: %d", nrow(X), ncol(X), length(unique(study))),
  "Same data, centering and study-level resampling as the registered method.", "",
  sprintf("Out-of-bag R2 : %.3f   (LASSO 0.040, Elastic Net 0.000)", r2_oob),
  sprintf("LOSO R2       : %.3f   <- the number that matters at n=%d", r2_loso, nrow(X)),
  sprintf("LOSO rho      : %+.3f", rho_loso),
  sprintf("Genera >= %.0f%% importance-stability: %d", 100*STABILITY_THR, sum(freq >= STABILITY_THR)),
  "", "Top genera by importance stability:")
for (i in seq_along(top)) out <- c(out, sprintf("  %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i]))
writeLines(out, file.path(R_OUT_DIR, "gdi", "random_forest.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "random_forest.txt"))
log_msg("=== DONE ===")
