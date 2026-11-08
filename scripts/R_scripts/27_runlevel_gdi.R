#!/usr/bin/env Rscript
# =============================================================================
# 27_runlevel_gdi.R
#
# A GDI built at RUN level from all 1,787 sequencing runs, with no GHI
# involved at any point. The target is the treatment label itself.
#
# RUN:  Rscript 27_runlevel_gdi.R
# TIME: 10-20 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS CAN USE EVERY RUN WHEN NOTHING ELSE COULD
# -----------------------------------------------------------------------------
# Every GHI-targeted analysis in this project collapses runs to arms, because
# a GHI value exists only at group level: thirty runs from one arm all share
# a single biomarker number, so treating them as thirty observations would be
# pseudo-replication - inflating apparent precision without adding
# information.
#
# The treatment LABEL does not have that problem. Every run from a challenged
# arm genuinely was a challenged bird. The label is a property of the run, not
# an average imposed on it. So the full 1,787 runs are legitimately available
# here, roughly twelve times the 141 arms every other model was limited to.
#
# What changes is the question. This does not ask "how healthy is this gut",
# it asks "was this bird challenged". Those are different questions, and a
# panel that answers the second is a challenge-detection classifier rather
# than a dysbiosis index in the clinical sense. It is arguably closer to the
# deployable product, since in the field the GHI is unknown too - but the
# distinction must be stated plainly rather than blurred.
#
# -----------------------------------------------------------------------------
# WHERE NON-INDEPENDENCE STILL BITES, AND HOW IT IS HANDLED
# -----------------------------------------------------------------------------
# Runs within a study are not independent - shared diet, house, extraction
# kit, sequencing batch. Ignoring that would inflate everything.
#
#   Stability selection resamples whole STUDIES at 50%, never individual runs.
#   Leave-one-study-out holds out every run of a study together.
#   Within-study centering removes each study's baseline before fitting.
#
# The honest estimate is therefore the LOSO AUC on held-out studies, not the
# in-sample fit. Both are reported so the gap is visible.
#
# -----------------------------------------------------------------------------
# WHICH RUNS ARE USED
# -----------------------------------------------------------------------------
# Control versus challenged only. Additive arms are excluded from the primary
# contrast because "not challenged" and "given a beneficial supplement" are
# different states, and pooling them would blur the thing being detected.
# Challenged-and-treated arms are excluded entirely - they are both.
# A secondary contrast including additives as a third class is reported after.
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
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")
HEALTHY   <- "Control"
CHALLENGE <- c("Pathogen_Challenge","Antibiotic_Challenge","Dietary_Feed_Challenge")

log_msg("[1/6] Loading the run-level CLR matrix")
clr <- fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv"))
M <- as.matrix(clr[, -1, with = FALSE]); rownames(M) <- clr[[names(clr)[1]]]
log_msg("  %d genera x %d runs", nrow(M), ncol(M))

md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), rid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), rid)

runs <- intersect(colnames(M), names(r2s))
M <- M[, runs, drop = FALSE]
log_msg("  %d runs matched to metadata", length(runs))

log_msg("[2/6] Attaching the treatment label (no GHI used)")
bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
bk <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||")
cls_lut <- tapply(trimws(as.character(bio$Challenge_Type)), bk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })
rk  <- paste(r2s[runs], r2t[runs], sep = "||")
cls <- unname(cls_lut[rk])
study_all <- unname(r2s[runs])

keep <- cls %in% c(HEALTHY, CHALLENGE)
Xr <- t(M[, keep, drop = FALSE])
yr <- as.integer(cls[keep] %in% CHALLENGE)   # 1 = challenged
sr <- study_all[keep]
log_msg("  %d runs retained: %d control, %d challenged, across %d studies",
        nrow(Xr), sum(yr == 0), sum(yr == 1), length(unique(sr)))

# a study must contain BOTH classes to contribute a within-study contrast
tb <- table(sr, yr)
both <- rownames(tb)[tb[, "0"] > 0 & tb[, "1"] > 0]
ok <- sr %in% both
Xr <- Xr[ok, , drop = FALSE]; yr <- yr[ok]; sr <- sr[ok]
log_msg("  %d studies contain both classes; %d runs retained", length(both), nrow(Xr))
log_msg("  (for comparison, every GHI-targeted model in this project used 141 ARMS)")

log_msg("[3/6] Within-study centering")
Xc <- Xr
for (s in unique(sr)) { i <- sr == s
  Xc[i, ] <- sweep(Xr[i, , drop = FALSE], 2, colMeans(Xr[i, , drop = FALSE]), "-") }
keep_g <- apply(Xc, 2, function(v) stats::sd(v) > 0)
Xc <- Xc[, keep_g, drop = FALSE]
log_msg("  %d genera retained after removing zero-variance columns", ncol(Xc))

log_msg("[4/6] Stability selection, %d resamples, %.0f%% of studies each", N_BOOT, 100*SUBFRAC)
studies <- unique(sr); ntake <- max(3, floor(length(studies) * SUBFRAC))
cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
for (b in seq_len(N_BOOT)) {
  tk <- sample(studies, ntake); idx <- which(sr %in% tk)
  if (length(idx) < 20 || length(unique(yr[idx])) < 2) next
  f <- tryCatch(cv.glmnet(Xc[idx, , drop = FALSE], yr[idx], family = "binomial",
                          alpha = 1, nfolds = 5, standardize = TRUE),
               error = function(e) NULL)
  if (is.null(f)) next
  co <- as.matrix(coef(f, s = f$lambda.1se)); co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  sel <- rownames(co)[co[, 1] != 0]
  cnt[sel] <- cnt[sel] + 1; nok <- nok + 1
}
freq <- if (nok > 0) cnt / nok else cnt
freq <- freq[!(names(freq) %in% UNASSIGNABLE)]
log_msg("  %d of %d resamples usable", nok, N_BOOT)
log_msg("  genera at the registered %.0f%% threshold: %d", 100*STABILITY_THR, sum(freq >= STABILITY_THR))
top <- sort(freq, decreasing = TRUE)[1:min(12, length(freq))]
for (i in seq_along(top)) log_msg("    %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i])

log_msg("[5/6] Leave-one-study-out - the honest test")
pred <- rep(NA_real_, length(yr))
for (s in unique(sr)) {
  tr <- sr != s; te <- sr == s
  if (length(unique(yr[tr])) < 2 || sum(te) < 2) next
  f <- tryCatch(cv.glmnet(Xc[tr, , drop = FALSE], yr[tr], family = "binomial",
                          alpha = 1, nfolds = 5, standardize = TRUE),
               error = function(e) NULL)
  if (is.null(f)) next
  pred[te] <- as.numeric(predict(f, newx = Xc[te, , drop = FALSE],
                                 s = f$lambda.1se, type = "response"))
}
okp <- !is.na(pred)
auc_of <- function(p, y) {
  p1 <- p[y == 1]; p0 <- p[y == 0]
  if (length(p1) < 3 || length(p0) < 3) return(NA_real_)
  r <- rank(c(p1, p0))
  (sum(r[seq_along(p1)]) - length(p1)*(length(p1)+1)/2) / (length(p1)*length(p0))
}
auc_loso <- auc_of(pred[okp], yr[okp])
log_msg("  LOSO AUC = %.3f  (on %d of %d runs, %d studies)",
        auc_loso, sum(okp), length(yr), length(unique(sr[okp])))
log_msg("  0.5 = chance; 0.70-0.80 is the usual minimum for a usable test.")

# in-sample, for the gap
f_full <- cv.glmnet(Xc, yr, family = "binomial", alpha = 1, nfolds = 5, standardize = TRUE)
p_in <- as.numeric(predict(f_full, newx = Xc, s = f_full$lambda.1se, type = "response"))
auc_in <- auc_of(p_in, yr)
n1se <- sum(as.matrix(coef(f_full, s = f_full$lambda.1se))[-1, 1] != 0)
log_msg("  in-sample AUC = %.3f, lambda.1se keeps %d genera", auc_in, n1se)
if (!is.na(auc_in) && !is.na(auc_loso) && auc_in - auc_loso > 0.15)
  log_msg("  *** large in-sample/LOSO gap: the model is fitting study-specific structure.")

log_msg("[6/6] Writing")
out <- c("RUN-LEVEL GDI - TREATMENT LABEL AS TARGET, NO GHI", strrep("=", 74), "",
  "Built from every sequencing run rather than arm means. Legitimate here",
  "because the treatment label is a property of each run, unlike a GHI value",
  "which exists only as a group average.", "",
  "This asks 'was this bird challenged', not 'how healthy is this gut'.",
  "A panel answering it is a challenge-detection classifier, not a dysbiosis",
  "index in the clinical sense.", "",
  sprintf("Runs: %d (%d control, %d challenged)", nrow(Xc), sum(yr==0), sum(yr==1)),
  sprintf("Studies contributing both classes: %d", length(both)),
  sprintf("Genera: %d", ncol(Xc)),
  sprintf("Genera at %.0f%% stability: %d", 100*STABILITY_THR, sum(freq >= STABILITY_THR)),
  sprintf("In-sample AUC: %.3f (%d genera at lambda.1se)", auc_in, n1se),
  sprintf("LOSO AUC: %.3f   <- the honest figure", auc_loso), "",
  "Strongest candidates:")
for (i in seq_along(top)) out <- c(out, sprintf("  %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i]))
writeLines(out, file.path(R_OUT_DIR, "gdi", "runlevel_gdi.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "runlevel_gdi.txt"))
log_msg("=== DONE ===")
