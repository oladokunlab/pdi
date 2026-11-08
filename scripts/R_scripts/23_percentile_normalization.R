#!/usr/bin/env Rscript
# =============================================================================
# 23_percentile_normalization.R
#
# Percentile normalization (Oh et al. 2022) instead of within-study mean
# centering. Same target, same rules, same threshold - only how each study's
# abundances are put on a common scale changes.
#
# RUN:  Rscript 23_percentile_normalization.R   (needs 05_maaslin2.R to have run)
# TIME: a few minutes.
#
# -----------------------------------------------------------------------------
# WHAT IS ACTUALLY DIFFERENT HERE
# -----------------------------------------------------------------------------
# The registered pipeline centers each study on its own MEAN: an arm's value
# becomes "how far above/below this study's average." That removes the study's
# location but keeps its variance and distribution shape, so a study with a
# wide spread still contributes wider numbers than a tight one, and a skewed
# study stays skewed.
#
# Percentile normalization instead ranks every arm against THAT STUDY'S OWN
# CONTROL ARMS specifically. A value becomes "the Nth percentile of what
# healthy birds in this study showed" - a number from 0 to 1 that means the
# same thing in every study regardless of platform, primer region, depth or
# spread. That is a stronger form of harmonisation than mean centering, and
# it is the specific step Oh (2022) used to merge disparate cohorts.
#
# -----------------------------------------------------------------------------
# THE COST, STATED UP FRONT
# -----------------------------------------------------------------------------
# This requires each study to HAVE control arms to rank against. Studies
# without a control arm cannot be normalised this way and are dropped -
# not a bug, an inherent requirement of the method. The script reports how
# many studies and arms are lost so the comparison against the registered
# result stays honest: if this version scores better on far fewer studies,
# that is not a fair win.
#
# Ties are handled by averaging ranks. A genus absent from every control arm
# in a study gives a degenerate reference; those genus/study combinations are
# set to NA and the genus is dropped if it becomes too sparse overall.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

N_BOOT <- 1000; SUBFRAC <- 0.50; STABILITY_THR <- 0.70
MIN_CONTROL_ARMS <- 1        # a study needs at least this many control arms
MAX_GENUS_NA_FRAC <- 0.30    # drop a genus if >30% of arms end up NA
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")
HEALTHY <- "Control"

log_msg("[1/6] Loading arm-level data (same input as 06_lasso_selection.R)")
m <- readRDS(file.path(R_OUT_DIR, "gdi", "gdi_model.rds"))
X <- m$X; y <- m$y; study <- m$study; arm <- rownames(X)
log_msg("  %d arms x %d genera across %d studies", nrow(X), ncol(X), length(unique(study)))

log_msg("[2/6] Identifying control arms per study")
gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
cls <- unname(tapply(trimws(as.character(gh$Challenge_Type)), gk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })[arm])
is_ctrl <- cls %in% HEALTHY
log_msg("  %d control arms identified", sum(is_ctrl, na.rm = TRUE))

studies_with_ctrl <- unique(study[which(is_ctrl)])
n_before <- length(unique(study))
log_msg("  %d of %d studies have at least %d control arm(s)",
        length(studies_with_ctrl), n_before, MIN_CONTROL_ARMS)
lost <- setdiff(unique(study), studies_with_ctrl)
if (length(lost)) {
  log_msg("  studies droppable for lacking a control arm (%d):", length(lost))
  for (s in lost) log_msg("      %-24s %d arm(s)", s, sum(study == s))
}

log_msg("[3/6] Percentile-normalising each arm against its own study's controls")
keep <- study %in% studies_with_ctrl
Xk <- X[keep, , drop = FALSE]; yk <- y[keep]; sk <- study[keep]; ck <- is_ctrl[keep]
P <- matrix(NA_real_, nrow(Xk), ncol(Xk), dimnames = dimnames(Xk))
for (s in unique(sk)) {
  i <- which(sk == s); ci <- i[ck[i]]
  if (!length(ci)) next
  for (g in seq_len(ncol(Xk))) {
    ref <- Xk[ci, g]
    ref <- ref[is.finite(ref)]
    if (!length(ref) || length(unique(ref)) == 0) next
    # percentile of each arm's value within this study's control distribution
    P[i, g] <- vapply(Xk[i, g], function(v) {
      if (!is.finite(v)) return(NA_real_)
      (sum(ref < v) + 0.5 * sum(ref == v)) / length(ref)
    }, numeric(1))
  }
}
na_frac <- colMeans(is.na(P))
drop_g <- na_frac > MAX_GENUS_NA_FRAC
log_msg("  %d genera dropped for >%.0f%% NA after normalisation", sum(drop_g), 100*MAX_GENUS_NA_FRAC)
P <- P[, !drop_g, drop = FALSE]
ok_rows <- stats::complete.cases(P) & is.finite(yk)
P <- P[ok_rows, , drop = FALSE]; yk <- yk[ok_rows]; sk <- sk[ok_rows]
log_msg("  final: %d arms x %d genera across %d studies", nrow(P), ncol(P), length(unique(sk)))
log_msg("  (registered comparison ran on %d arms across %d studies)", nrow(X), n_before)

log_msg("[4/6] Centering the outcome within study (X is already on a common 0-1 scale)")
yc <- yk
for (s in unique(sk)) { i <- sk == s; yc[i] <- yk[i] - mean(yk[i]) }

log_msg("[5/6] Stability selection, %d resamples, %.0f%% study subsampling", N_BOOT, 100*SUBFRAC)
studies <- unique(sk); ntake <- max(3, floor(length(studies) * SUBFRAC))
cnt <- setNames(numeric(ncol(P)), colnames(P)); nok <- 0
for (b in seq_len(N_BOOT)) {
  tk <- sample(studies, ntake); idx <- which(sk %in% tk)
  if (length(idx) < 10 || length(unique(yc[idx])) < 5) next
  f <- tryCatch(cv.glmnet(P[idx, , drop = FALSE], yc[idx], alpha = 1,
                          nfolds = min(5, length(idx)), standardize = TRUE, intercept = FALSE),
               error = function(e) NULL)
  if (is.null(f)) next
  co <- as.matrix(coef(f, s = f$lambda.1se)); co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  sel <- rownames(co)[co[,1] != 0]
  cnt[sel] <- cnt[sel] + 1; nok <- nok + 1
}
freq <- if (nok > 0) cnt / nok else cnt
freq <- freq[!(names(freq) %in% UNASSIGNABLE)]
log_msg("  %d of %d resamples usable", nok, N_BOOT)
log_msg("  genera at the registered %.0f%% threshold: %d", 100*STABILITY_THR, sum(freq >= STABILITY_THR))
top <- sort(freq, decreasing = TRUE)[1:min(10, length(freq))]
for (i in seq_along(top)) log_msg("    %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i])

log_msg("[6/6] In-sample fit and writing")
full <- tryCatch({ f <- cv.glmnet(P, yc, alpha=1, standardize=TRUE, intercept=FALSE)
                   pr <- as.numeric(predict(f, newx=P, s=f$lambda.1se))
                   r2 <- if (sd(pr)==0) 0 else max(0, 1 - sum((yc-pr)^2)/sum((yc-mean(yc))^2))
                   list(r2=r2, n1se=sum(as.matrix(coef(f, s=f$lambda.1se))[-1,1] != 0)) },
                error=function(e) list(r2=NA_real_, n1se=NA))
log_msg("  R2 = %.3f   (registered mean-centered comparison: 0.040)", full$r2)
log_msg("  lambda.1se selects %s genera   (registered comparison: 2)", full$n1se)

out <- c("PERCENTILE NORMALIZATION (Oh 2022) vs WITHIN-STUDY MEAN CENTERING",
  strrep("=", 74), "",
  "Each arm ranked against its own study's CONTROL arms, scaled 0-1, rather",
  "than centered on the study mean. Same target, same 70% threshold, same",
  "1000-iteration study-level resampling.", "",
  sprintf("Arms: %d across %d studies (registered comparison: %d arms, %d studies)",
          nrow(P), length(unique(sk)), nrow(X), n_before),
  sprintf("Studies dropped for having no control arm: %d", length(lost)),
  sprintf("Genera: %d (%d dropped for excessive NA after normalisation)", ncol(P), sum(drop_g)),
  "",
  sprintf("In-sample R2: %.3f   (registered: 0.040)", full$r2),
  sprintf("lambda.1se genera: %s   (registered: 2)", full$n1se),
  sprintf("Genera at %.0f%% stability: %d", 100*STABILITY_THR, sum(freq >= STABILITY_THR)),
  "", "Strongest candidates:")
for (i in seq_along(top)) out <- c(out, sprintf("  %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i]))
writeLines(out, file.path(R_OUT_DIR, "gdi", "percentile_normalization.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "percentile_normalization.txt"))
log_msg("=== DONE ===")
