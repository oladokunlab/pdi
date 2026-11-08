#!/usr/bin/env Rscript
# =============================================================================
# 24_elastic_net.R
#
# Elastic Net (alpha < 1) instead of LASSO (alpha = 1). Everything else
# identical to 06_lasso_selection.R: same arms, same target, same
# within-study centering, same 1000-iteration study-level resampling, same
# registered 70% stability threshold.
#
# RUN:  Rscript 24_elastic_net.R    (needs 06_lasso_selection.R to have run)
# TIME: ~5 minutes (sweeps several alpha values).
#
# -----------------------------------------------------------------------------
# WHY THIS IS NOT JUST ANOTHER PARAMETER TWEAK
# -----------------------------------------------------------------------------
# LASSO's L1 penalty forces correlated predictors to COMPETE: given a group of
# genera that co-occur, it keeps one and zeros the rest. Which one it keeps
# depends on which studies happen to be in a given resample, so a genuinely
# real but correlated group splits its selection frequency across members and
# no single member reaches a high stability score.
#
# That is a plausible explanation for the pattern seen across every
# configuration tested in this project: max stability plateauing at 27-34%
# regardless of arms, target, or scaling. A ceiling that stable across
# unrelated changes looks more like a property of the estimator than of the
# data.
#
# Elastic Net's mixed L1/L2 penalty lets correlated predictors ENTER TOGETHER
# (the "grouping effect", Zou & Hastie 2005). If the ceiling is caused by
# vote-splitting, stability should rise as alpha decreases. If the ceiling is
# caused by there being little signal, it will not move - and that is an
# informative answer too, since it rules out the estimator as the explanation.
#
# Wirbel et al. (SIAMCAT, Genome Biology 2021) benchmarked 7,488 model
# configurations across 130 microbiome classification tasks and found Elastic
# Net outperformed LASSO on average (p = 0.001), which is the empirical basis
# for trying it here.
#
# -----------------------------------------------------------------------------
# WHAT IS REPORTED
# -----------------------------------------------------------------------------
# The full alpha sweep, so the trend is visible rather than a single number:
#   alpha = 1.00  pure LASSO (the registered method, for direct comparison)
#   alpha = 0.75, 0.50, 0.25  progressively more L2 / more grouping
#   alpha = 0.05  nearly ridge; keeps almost everything, useful as a bound
#
# Reported per alpha: max stability, how many genera clear 70%, in-sample R2,
# and the top genera. A rising stability curve as alpha falls is the
# signature of vote-splitting having been the constraint.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

ALPHAS <- c(1.00, 0.75, 0.50, 0.25, 0.05)
N_BOOT <- 1000; SUBFRAC <- 0.50; STABILITY_THR <- 0.70
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")

log_msg("[1/4] Loading the same arm-level data 06_lasso_selection.R used")
m <- readRDS(file.path(R_OUT_DIR, "gdi", "gdi_model.rds"))
X <- m$X; y <- m$y; study <- m$study
log_msg("  %d arms x %d genera across %d studies", nrow(X), ncol(X), length(unique(study)))

log_msg("[2/4] Within-study centering (identical to the registered method)")
Xc <- X; yc <- y
for (s in unique(study)) { i <- study == s
  Xc[i, ] <- sweep(X[i, , drop = FALSE], 2, colMeans(X[i, , drop = FALSE]), "-")
  yc[i] <- y[i] - mean(y[i]) }

log_msg("[3/4] Alpha sweep, %d resamples each", N_BOOT)
studies <- unique(study); ntake <- max(3, floor(length(studies) * SUBFRAC))
results <- list()

for (a in ALPHAS) {
  cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
  for (b in seq_len(N_BOOT)) {
    tk <- sample(studies, ntake); idx <- which(study %in% tk)
    if (length(idx) < 10 || length(unique(yc[idx])) < 5) next
    f <- tryCatch(cv.glmnet(Xc[idx, , drop = FALSE], yc[idx], alpha = a,
                            nfolds = min(5, length(idx)), standardize = TRUE, intercept = FALSE),
                 error = function(e) NULL)
    if (is.null(f)) next
    co <- as.matrix(coef(f, s = f$lambda.1se)); co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
    sel <- rownames(co)[co[,1] != 0]
    cnt[sel] <- cnt[sel] + 1; nok <- nok + 1
  }
  freq <- if (nok > 0) cnt / nok else cnt
  freq <- freq[!(names(freq) %in% UNASSIGNABLE)]

  fit <- tryCatch({ f <- cv.glmnet(Xc, yc, alpha = a, standardize = TRUE, intercept = FALSE)
                    pr <- as.numeric(predict(f, newx = Xc, s = f$lambda.1se))
                    r2 <- if (sd(pr) == 0) 0 else max(0, 1 - sum((yc-pr)^2)/sum((yc-mean(yc))^2))
                    n1se <- sum(as.matrix(coef(f, s = f$lambda.1se))[-1,1] != 0)
                    list(r2 = r2, n1se = n1se) },
                  error = function(e) list(r2 = NA_real_, n1se = NA))

  results[[as.character(a)]] <- list(alpha = a, freq = freq, r2 = fit$r2, n1se = fit$n1se,
                                     max_stab = if (length(freq)) max(freq) else NA,
                                     n70 = sum(freq >= STABILITY_THR))
  log_msg("  alpha=%.2f  max stability=%3.0f%%  >=70%%: %d  R2=%.3f  lambda.1se genera=%s",
          a, 100*results[[as.character(a)]]$max_stab, results[[as.character(a)]]$n70,
          fit$r2, fit$n1se)
}

log_msg("")
log_msg("[4/4] Top genera per alpha")
for (a in ALPHAS) {
  r <- results[[as.character(a)]]
  top <- sort(r$freq, decreasing = TRUE)[1:min(6, length(r$freq))]
  log_msg("  alpha=%.2f: %s", a,
          paste(sprintf("%s=%.0f%%", gsub("^g__","",names(top)), 100*top), collapse="  "))
}

best <- results[[which.max(vapply(results, function(z) ifelse(is.na(z$max_stab), 0, z$max_stab), numeric(1)))]]
log_msg("")
log_msg("Highest stability at alpha=%.2f (%.0f%%), versus %.0f%% for pure LASSO.",
        best$alpha, 100*best$max_stab, 100*results[["1"]]$max_stab)
if (best$max_stab > results[["1"]]$max_stab + 0.10) {
  log_msg("  A rise of this size suggests correlated genera were splitting their")
  log_msg("  selection frequency under L1 - the estimator, not the data, was the")
  log_msg("  limiting factor. Whether it clears the registered 70%% bar is separate.")
} else {
  log_msg("  Stability did not rise materially as alpha fell. Vote-splitting among")
  log_msg("  correlated predictors is therefore NOT the explanation for the ceiling;")
  log_msg("  this points back to limited signal rather than estimator choice.")
}

out <- c("ELASTIC NET ALPHA SWEEP vs LASSO", strrep("=", 74), "",
  sprintf("Arms: %d, genera: %d, studies: %d", nrow(X), ncol(X), length(unique(study))),
  "Identical data, centering, resampling and threshold to 06_lasso_selection.R.",
  "Only the glmnet alpha parameter changes (1.00 = the registered LASSO).", "",
  sprintf("  %-8s %12s %8s %8s %10s", "alpha", "max stab", ">=70%", "R2", "1se genera"))
for (a in ALPHAS) { r <- results[[as.character(a)]]
  out <- c(out, sprintf("  %-8.2f %11.0f%% %8d %8.3f %10s", a, 100*r$max_stab, r$n70, r$r2, r$n1se)) }
out <- c(out, "", "Top genera per alpha:")
for (a in ALPHAS) { r <- results[[as.character(a)]]
  top <- sort(r$freq, decreasing = TRUE)[1:min(6, length(r$freq))]
  out <- c(out, sprintf("  alpha=%.2f: %s", a,
           paste(sprintf("%s=%.0f%%", gsub("^g__","",names(top)), 100*top), collapse="  "))) }
writeLines(out, file.path(R_OUT_DIR, "gdi", "elastic_net_sweep.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "elastic_net_sweep.txt"))
log_msg("=== DONE ===")
