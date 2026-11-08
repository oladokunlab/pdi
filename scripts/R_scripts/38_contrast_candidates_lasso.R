#!/usr/bin/env Rscript
# =============================================================================
# 38_contrast_candidates_lasso.R
#
# Does the expanded candidate set from 37_per_taxon_contrast.R combine into a
# panel? Same LASSO + stability selection procedure as 06_lasso_selection.R,
# restricted to genera that showed a significant contrast in 37's output.
#
# RUN:  Rscript 38_contrast_candidates_lasso.R
#
# -----------------------------------------------------------------------------
# WHY THIS IS THE RIGHT NEXT TEST
# -----------------------------------------------------------------------------
# 37_per_taxon_contrast.R tests each of 84 genera ONE AT A TIME against the
# GHI, using whichever studies have a control arm and a treated arm. Several
# genera passed p <= 0.05 in at least one contrast, with far more studies
# contributing than the correlation-based test in 18_taxon_metaanalysis.R
# could use. That answers "is there a candidate association" - it does not
# answer "do these combine into something that predicts health on a study
# the model has never seen". Only stability selection under LOSO answers
# that, which is what 06_lasso_selection.R already tested on the full
# 84-genus set and found nothing at the 70% threshold.
#
# This script asks the same question again, but starting from a smaller,
# pre-screened set rather than all 84, to see whether restricting to genera
# that already show SOME signal changes the answer.
#
# -----------------------------------------------------------------------------
# A CAVEAT THAT MUST BE STATED PLAINLY, NOT BURIED
# -----------------------------------------------------------------------------
# The candidate genera were selected because they correlate with the GHI in
# 37's contrast test. This script then tests whether they predict the SAME
# GHI in LASSO. That is not two independent pieces of evidence - it is the
# same outcome variable used twice, once to choose the candidates and once
# to test them. This is a known form of selection bias (sometimes called
# double dipping or circular analysis): pre-filtering on significance before
# a second test on the same outcome makes finding "significance" easier than
# it should be, because chance associations in 37 get a second, biased shot
# at surviving here.
#
# The registered method deliberately avoided this by testing all 84 genera
# in LASSO without pre-filtering by univariate significance against the same
# target. That remains the methodologically correct approach and the result
# already reported: 0 genera at 70%, R2 = 0.000.
#
# If THIS script finds a panel where the registered method did not, that
# finding should be treated as a demonstration of what pre-filtering makes
# look possible, not as a more reliable discovery, and should not be reported
# as validated without a fully independent replication.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
suppressPackageStartupMessages(library(glmnet))

N_BOOT <- 1000
SUBSAMPLE_FRAC <- 0.50
STABILITY_THR <- 0.70
set.seed(42)

log_msg("[1/5] Loading candidate genera from 37's output")
contrast_tab <- fread(file.path(R_OUT_DIR, "gdi", "taxon_contrast_meta.tsv"))
candidates <- unique(contrast_tab$genus[!is.na(contrast_tab$p) & contrast_tab$p <= 0.05])
log_msg("  %d genera pass p <= 0.05 in at least one contrast (ADDITIVE/CHALLENGE/COMBINED)", length(candidates))
for (g in sort(candidates)) log_msg("      %s", g)

log_msg("")
log_msg("[2/5] Rebuilding the centered arm-level matrix directly from source data")
log_msg("  NOTE: gdi_model.rds was checked first and found to hold only 22 genus")
log_msg("  columns, not the 84 the registered LASSO is documented as using. Rather")
log_msg("  than guess which is right, the matrix is rebuilt here from")
log_msg("  genus_clr_batch_adjusted.tsv directly, the same primary source every")
log_msg("  other script in this project reads from.")

clr <- as.data.frame(fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")))
gname_col <- names(clr)[1]
Xr <- as.matrix(clr[, -1, drop = FALSE]); rownames(Xr) <- clr[[gname_col]]
runs_all <- colnames(Xr)

md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2study <- setNames(trimws(as.character(md$Study_ID)), rid)
r2arm   <- setNames(trimws(as.character(md$Treatment_Group)), rid)

ghi <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
arm_ghi <- setNames(ghi$GHI_Combined_raw,
                     paste(trimws(as.character(ghi$Study_ID)),
                           trimws(as.character(ghi$Treatment_Group)), sep = "||"))

common <- intersect(runs_all, names(r2study))
arm_key <- paste(r2study[common], r2arm[common], sep = "||")
uniq_arms <- unique(arm_key)
Aall <- t(vapply(uniq_arms, function(a) rowMeans(Xr[, common[arm_key == a], drop = FALSE], na.rm = TRUE),
                  numeric(nrow(Xr))))
colnames(Aall) <- rownames(Xr)
arm_study_all <- sub("\\|\\|.*$", "", uniq_arms)
y_all <- unname(arm_ghi[uniq_arms])

ok <- !is.na(y_all)
Aall <- Aall[ok, , drop = FALSE]; arm_study_all <- arm_study_all[ok]; y_all <- y_all[ok]

# drop single-arm studies (cannot contribute to within-study centering)
tab_s <- table(arm_study_all)
keep_s <- arm_study_all %in% names(tab_s)[tab_s >= 2]
Aall <- Aall[keep_s, , drop = FALSE]; arm_study_all <- arm_study_all[keep_s]; y_all <- y_all[keep_s]

# within-study centering, both X and y
X_full <- Aall
for (s in unique(arm_study_all)) {
  idx <- arm_study_all == s
  X_full[idx, ] <- sweep(Aall[idx, , drop = FALSE], 2, colMeans(Aall[idx, , drop = FALSE]), "-")
  y_all[idx] <- y_all[idx] - mean(y_all[idx])
}
y <- y_all
study <- arm_study_all
log_msg("  full matrix: %d arms x %d genera (all genera present, all 41 studies)",
        nrow(X_full), ncol(X_full))

keep <- intersect(candidates, colnames(X_full))
if (length(keep) < length(candidates)) {
  log_msg("  NOTE: %d candidate genus name(s) did not match and were dropped: %s",
          length(candidates) - length(keep),
          paste(setdiff(candidates, keep), collapse = ", "))
}
X <- X_full[, keep, drop = FALSE]
log_msg("  restricted matrix: %d arms x %d candidate genera", nrow(X), ncol(X))
if (ncol(X) < 2) {
  log_msg("FATAL: fewer than 2 candidate genera matched - cannot fit LASSO. Stopping.")
  quit(status = 1)
}

log_msg("")
log_msg("[3/5] LASSO with stability selection, identical procedure to 06_lasso_selection.R")
log_msg("  N_BOOT=%d, SUBSAMPLE_FRAC=%.2f, STABILITY_THR=%.2f", N_BOOT, SUBSAMPLE_FRAC, STABILITY_THR)

studies_u <- unique(study)
sel_count <- setNames(numeric(ncol(X)), colnames(X))
n_ok <- 0
r2_seeds <- numeric(0)

for (b in seq_len(N_BOOT)) {
  s_sub <- sample(studies_u, size = round(length(studies_u) * SUBSAMPLE_FRAC))
  idx <- study %in% s_sub
  Xb <- X[idx, , drop = FALSE]; yb <- y[idx]; sb <- study[idx]
  if (length(unique(sb)) < 3 || sum(idx) < 6) next
  fit <- tryCatch(
    cv.glmnet(Xb, yb, alpha = 1, nfolds = min(5, length(unique(sb))),
              intercept = FALSE, standardize = TRUE),
    error = function(e) NULL)
  if (is.null(fit)) next
  cf <- as.matrix(coef(fit, s = "lambda.1se"))[-1, 1]
  sel_count[names(cf)[cf != 0]] <- sel_count[names(cf)[cf != 0]] + 1
  n_ok <- n_ok + 1
}
stability <- sel_count / n_ok
log_msg("  %d of %d resamples usable", n_ok, N_BOOT)

log_msg("")
log_msg("[4/5] R2, median across 20 CV seeds (same diagnostic as the registered run)")
for (sd in 1:20) {
  set.seed(sd)
  fit <- tryCatch(cv.glmnet(X, y, alpha = 1, intercept = FALSE, standardize = TRUE),
                   error = function(e) NULL)
  if (is.null(fit)) { r2_seeds <- c(r2_seeds, NA); next }
  pred <- predict(fit, newx = X, s = "lambda.1se")
  ss_res <- sum((y - pred)^2); ss_tot <- sum((y - mean(y))^2)
  r2_seeds <- c(r2_seeds, max(0, 1 - ss_res / ss_tot))
}
log_msg("  R2 median = %.3f (range %.3f-%.3f across seeds)",
        median(r2_seeds, na.rm=TRUE), min(r2_seeds, na.rm=TRUE), max(r2_seeds, na.rm=TRUE))

log_msg("")
log_msg("=========================================================")
log_msg("[5/5] STABILITY RESULTS")
log_msg("=========================================================")
ord <- order(-stability)
log_msg("  %-28s %10s", "genus", "stability")
for (g in names(stability)[ord]) {
  st <- stability[g]
  tag <- if (!is.na(st) && st >= STABILITY_THR) "  *** PASSES 70% ***" else ""
  log_msg("  %-28s %9s%s", g, if (is.na(st)) "NA" else sprintf("%.1f%%", st*100), tag)
}

n_pass <- sum(stability >= STABILITY_THR, na.rm = TRUE)
log_msg("")
log_msg("  Genera at >= 70%% stability: %d", n_pass)
log_msg("  Max stability achieved: %.1f%% (%s)", max(stability)*100, names(which.max(stability)))
log_msg("  R2: %.3f", median(r2_seeds, na.rm=TRUE))

out <- data.frame(genus = names(stability), stability = stability)
out <- out[order(-out$stability), ]
fwrite(out, file.path(R_OUT_DIR, "gdi", "contrast_candidates_lasso_stability.tsv"), sep = "\t")

log_msg("")
log_msg("=== DONE ===")
log_msg("")
log_msg("READING: this result carries the double-dipping caveat stated at the top of")
log_msg("this script. The candidate genera were chosen because they already correlate")
log_msg("with the GHI in 37's output; testing them again against the same GHI here is")
log_msg("not independent confirmation. Compare against the registered result (all 84")
log_msg("genera, no pre-filtering): 0 at 70%%, R2 = 0.000.")
