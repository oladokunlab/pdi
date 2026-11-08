#!/usr/bin/env Rscript
# =============================================================================
# 08_loso_validate.R
#
# Leave-one-study-out validation of the GDI. The only number in this project
# that estimates performance on a study the model has never seen.
#
# RUN:  Rscript 08_loso_validate.R      (needs 05, 06, 07 first)
# TIME: 15-45 minutes.
#
# -----------------------------------------------------------------------------
# WHY THE PANEL IS RE-SELECTED INSIDE EVERY FOLD
# -----------------------------------------------------------------------------
# The panel in 07 was selected using ALL studies. Holding out a study but
# scoring it with that panel asks "how well does a model built partly from this
# study predict this study" - the held-out data influenced which genera were
# chosen. That is selection leakage and it inflates cross-study performance.
#
# Here LASSO selection and centroid construction are re-run from scratch inside
# each fold using only the training studies. The held-out study contributes
# nothing to selection, nothing to the centroids, nothing to the median split.
#
# RESIDUAL LEAKAGE, DISCLOSED: the MaAsLin2 pre-filter in 05 ran once on the
# full corpus, so the candidate pool entering each fold was informed by all
# studies. Re-running it per fold costs hours. The pre-filter narrows ~84
# genera to ~29, so it constrains rather than determines the panel, but the
# reported performance is slightly optimistic and must be described that way.
#
# -----------------------------------------------------------------------------
# TWO SCORING SPACES, AND WHY BOTH MATTER
# -----------------------------------------------------------------------------
# CENTERED  - the held-out study's arms are centred on that study's own
#             abundance means before scoring. Matches the space the panel was
#             selected in. Answers "does the panel carry real signal". Requires
#             several samples from the same study, so it describes a
#             within-flock comparison.
#
# ABSOLUTE  - raw CLR values, no centring. One sample can be scored with no
#             reference group, which is what a deployed qPCR assay must do.
#             Exposed to between-study technical variation and expected to be
#             weaker: in 07 it gave r = -0.235 against -0.455 centred.
#
# CENTERED is the scientific claim; ABSOLUTE is the deployment claim and is the
# honest number to quote for a qPCR panel. Panel SELECTION stays in centred
# space for both, because the uncentred fit was diagnosed in 06 as having no
# signal to select from (CV error 151.96 vs null variance 151.49).
#
# Centring the held-out study uses its ABUNDANCES only, never its GHI, so no
# outcome information crosses the fold boundary.
#
# -----------------------------------------------------------------------------
# WHY A PERMUTATION TEST AND NOT A SIGN TEST
# -----------------------------------------------------------------------------
# Studies contribute 2-7 arms. A per-study Spearman on 3 arms takes only a few
# possible values and is frequently exactly +/-1, so counting how many studies
# point the right way is close to counting coin flips - the earlier sign test
# (18/29, p = 0.13) was underpowered by construction, not because the index
# failed.
#
# The pooled correlation uses every arm, but its ordinary p-value is invalid
# because arms within a study are not independent. The test below permutes the
# outcome WITHIN each study and recomputes the pooled statistic, preserving the
# clustering and giving a valid p-value for the statistic actually reported.
# Per-study values are still printed, for inspection rather than inference.
#
# -----------------------------------------------------------------------------
# EXPECTED PERFORMANCE
# -----------------------------------------------------------------------------
# Continuous outcome, so correlation is the metric (Asnicar 2025: Spearman
# 0.30-0.46 across independent validation cohorts for a comparable design).
# AUC against the curated label is an external check - that label is never used
# in construction. Cross-study drops of 5-15 points versus within-study fits
# are normal; Pietruska & Hauck saw 0.71-0.99 fall to 0.42-0.68.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

N_BOOT_FOLD     <- 200
STABILITY_FLOOR <- 0.20
SUBSAMPLE_FRAC  <- 0.50
LAMBDA_RULE     <- { e <- Sys.getenv("LAMBDA_RULE", ""); if (nzchar(e)) e else "1se" }
MIN_PANEL       <- 3
N_PERM          <- 5000

# Placeholder taxonomy: not organisms, cannot be qPCR targets. Must match 07.
UNASSIGNABLE <- c("g__uncultured", "g__unclassified", "g__uncultured_bacterium",
                  "g__metagenome", "g__Incertae_Sedis", "g__unidentified")

# Contrast definitions, identical to 04e so both validations speak about the
# same groups.
HEALTHY   <- c("Control")
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")
ADDITIVE  <- c("Non_Disruptive_Additive")

set.seed(42)

log_msg("[1/7] Loading")
mf <- file.path(R_OUT_DIR, "gdi", "gdi_model.rds")
if (!file.exists(mf)) stop("Not found: ", mf, "\n  Run 06_lasso_selection.R first.")
m <- readRDS(mf)
X_all <- m$X; y_all <- m$y; study_all <- m$study
if (is.null(X_all) || is.null(y_all)) stop("gdi_model.rds is missing X/y - re-run 06.")
arm_names <- rownames(X_all)
log_msg("  %d arms, %d genera, %d studies", nrow(X_all), ncol(X_all),
        length(unique(study_all)))

# ---- curated contrast label -------------------------------------------------
# Read Challenge_Type from the GHI score table rather than guessing from
# treatment-group strings. An earlier version pattern-matched labels like "CK",
# "BL" and "0 mg/kg", could classify only 40 of 135 arms, and any mistake there
# silently inverts the AUC. The curated column is authoritative.
arm_class <- rep(NA_character_, length(arm_names))
gtsv <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
if (file.exists(gtsv)) {
  gh <- fread(gtsv)
  if ("Challenge_Type" %in% names(gh)) {
    k_gh <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
    lut <- tapply(trimws(as.character(gh$Challenge_Type)), k_gh, function(v) {
      v <- v[!is.na(v) & nzchar(v)]
      if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1]
    })
    arm_class <- unname(lut[arm_names])
    log_msg("  contrast label from Challenge_Type: %d of %d arms",
            sum(!is.na(arm_class)), length(arm_names))
  } else log_msg("  NOTE: Challenge_Type absent from %s - AUC skipped", basename(gtsv))
} else log_msg("  NOTE: %s not found - AUC skipped", basename(gtsv))

log_msg("  control=%d  challenged=%d  additive=%d  other/unlabelled=%d",
        sum(arm_class %in% HEALTHY), sum(arm_class %in% CHALLENGE),
        sum(arm_class %in% ADDITIVE), sum(is.na(arm_class)))

dist_to <- function(row, cen) sqrt(sum((row - cen)^2))
get_coefs <- function(fit, l) {
  co <- as.matrix(coef(fit, s = l))
  co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  co[co[, 1] != 0, 1]
}

# =============================================================================
# STEP 2 - FOLD LOOP
# =============================================================================
log_msg("[2/7] Leave-one-study-out (%d folds; selection centred, scoring both spaces)",
        length(unique(study_all)))

studies <- sort(unique(study_all))
rows <- list(); fold_panels <- list(); n_fallback <- 0

for (s_out in studies) {
  tr <- study_all != s_out; te <- !tr
  if (sum(te) < 2) { log_msg("  %-24s skipped (1 arm)", s_out); next }

  Xtr <- X_all[tr, , drop = FALSE]; ytr <- y_all[tr]; str_ <- study_all[tr]
  keep <- str_ %in% names(table(str_))[table(str_) >= 2]
  Xtr <- Xtr[keep, , drop = FALSE]; ytr <- ytr[keep]; str_ <- str_[keep]
  if (length(unique(str_)) < 5) { log_msg("  %-24s skipped (training too small)", s_out); next }

  Xtrc <- Xtr; ytrc <- ytr
  for (s in unique(str_)) {
    i <- str_ == s
    Xtrc[i, ] <- sweep(Xtr[i, , drop = FALSE], 2, colMeans(Xtr[i, , drop = FALSE]), "-")
    ytrc[i]   <- ytr[i] - mean(ytr[i])
  }

  st_tr <- unique(str_)
  n_take <- max(3, floor(length(st_tr) * SUBSAMPLE_FRAC))
  counts <- setNames(numeric(ncol(Xtrc)), colnames(Xtrc)); n_ok <- 0
  for (b in seq_len(N_BOOT_FOLD)) {
    take <- sample(st_tr, n_take); idx <- which(str_ %in% take)
    if (length(idx) < 10 || length(unique(ytrc[idx])) < 5) next
    fit <- tryCatch(cv.glmnet(Xtrc[idx, , drop = FALSE], ytrc[idx], alpha = 1,
                              nfolds = min(5, length(idx)), standardize = TRUE,
                              intercept = FALSE), error = function(e) NULL)
    if (is.null(fit)) next
    sel <- get_coefs(fit, if (LAMBDA_RULE == "1se") fit$lambda.1se else fit$lambda.min)
    counts[names(sel)] <- counts[names(sel)] + 1; n_ok <- n_ok + 1
  }
  if (n_ok < 20) { log_msg("  %-24s skipped (%d usable resamples)", s_out, n_ok); next }
  freq <- counts / n_ok

  full <- tryCatch({
    cvf <- cv.glmnet(Xtrc, ytrc, alpha = 1, standardize = TRUE, intercept = FALSE)
    get_coefs(cvf, if (LAMBDA_RULE == "1se") cvf$lambda.1se else cvf$lambda.min)
  }, error = function(e) numeric(0))

  panel_f <- setdiff(names(freq)[freq >= STABILITY_FLOOR & names(freq) %in% names(full)],
                     UNASSIGNABLE)
  fallback <- FALSE
  if (length(panel_f) < MIN_PANEL) {
    panel_f <- setdiff(names(sort(freq, decreasing = TRUE)), UNASSIGNABLE)[seq_len(MIN_PANEL)]
    fallback <- TRUE; n_fallback <- n_fallback + 1
  }
  fold_panels[[s_out]] <- panel_f

  yte <- y_all[te]
  rec <- data.table(study = s_out, arm = arm_names[te], ghi = yte,
                    ghi_centered = yte - mean(yte), cls = arm_class[te],
                    n_panel = length(panel_f), fallback = fallback,
                    gdi_centered = NA_real_, gdi_absolute = NA_real_)

  for (sp in c("centered", "absolute")) {
    Xp <- if (sp == "centered") Xtrc[, panel_f, drop = FALSE] else Xtr[, panel_f, drop = FALSE]
    yv <- if (sp == "centered") ytrc else ytr
    hh <- yv > stats::median(yv)
    if (sum(hh) < 3 || sum(!hh) < 3) next
    ch <- colMeans(Xp[hh, , drop = FALSE]); cd <- colMeans(Xp[!hh, , drop = FALSE])
    Xte <- X_all[te, panel_f, drop = FALSE]
    if (sp == "centered") Xte <- sweep(Xte, 2, colMeans(Xte), "-")
    # HIGH = MORE DYSBIOTIC, matching 07. Keep these two identical.
    set(rec, j = paste0("gdi_", sp),
        value = apply(Xte, 1, function(r) dist_to(r, ch) - dist_to(r, cd)))
  }

  rows[[s_out]] <- rec
  rr <- if (sum(!is.na(rec$gdi_centered)) >= 3 &&
            uniqueN(rec$gdi_centered) > 1 && uniqueN(rec$ghi_centered) > 1)
          suppressWarnings(cor(rec$gdi_centered, rec$ghi_centered, method = "spearman")) else NA_real_
  log_msg("  %-24s arms=%2d panel=%2d%s rho=%s", s_out, sum(te), length(panel_f),
          if (fallback) " (fb)" else "     ",
          if (is.na(rr)) "  n/a" else sprintf("%+.3f", rr))
}

if (!length(rows)) stop("No fold produced a prediction.")
pred <- rbindlist(rows, fill = TRUE)

# =============================================================================
# STEP 3 - POOLED CORRELATION WITH A WITHIN-STUDY PERMUTATION TEST
# =============================================================================
log_msg("[3/7] Pooled correlation (negative = working direction)")
log_msg("  %d arms, %d studies, %d fold(s) used a fallback panel",
        nrow(pred), uniqueN(pred$study), n_fallback)

perm_p <- function(g, y, grp, n = N_PERM) {
  ok <- !is.na(g) & !is.na(y)
  g <- g[ok]; y <- y[ok]; grp <- grp[ok]
  if (length(g) < 10) return(list(rho = NA_real_, p = NA_real_))
  obs <- suppressWarnings(cor(g, y, method = "spearman"))
  idx_by <- split(seq_along(g), grp)
  cnt <- 0L
  for (i in seq_len(n)) {
    yp <- y
    for (j in idx_by) if (length(j) > 1) yp[j] <- y[sample(j)]
    if (suppressWarnings(cor(g, yp, method = "spearman")) <= obs) cnt <- cnt + 1L
  }
  list(rho = obs, p = (cnt + 1) / (n + 1))
}

pool <- list()
# Each space is tested against the outcome that MATCHES it. Comparing an
# uncentred GDI against a within-study-centred GHI mixes spaces: the predictor
# retains between-study variation that the outcome has had removed by
# construction, so the correlation is diluted for a mechanical reason rather
# than a biological one. The centred index is judged on within-study contrasts;
# the absolute index is judged on the raw GHI, which is what a deployed assay
# would actually have to predict.
OUTCOME <- list(centered = "ghi_centered", absolute = "ghi")
for (sp in c("centered", "absolute")) {
  r <- perm_p(pred[[paste0("gdi_", sp)]], pred[[OUTCOME[[sp]]]], pred$study)
  pool[[sp]] <- r
  log_msg("  %-9s vs %-13s rho = %+.3f   permutation p = %.4f",
          sp, OUTCOME[[sp]], r$rho, r$p)
}
# Also reported: absolute scored against the centred outcome, which is the
# comparison an earlier version made. Kept so the two are not confused.
r_x <- perm_p(pred$gdi_absolute, pred$ghi_centered, pred$study)
log_msg("  absolute  vs ghi_centered  rho = %+.3f  p = %.4f  (mismatched spaces,",
        r_x$rho, r_x$p)
log_msg("            reported only to show it is not the deployment number)")
log_msg("  GHI is shuffled WITHIN each study, preserving the clustering that")
log_msg("  makes an ordinary correlation p-value invalid here.")

# =============================================================================
# STEP 4 - PER-STUDY (INSPECTION, NOT INFERENCE)
# =============================================================================
log_msg("[4/7] Per-study Spearman, centred space")
per <- pred[, .(n = .N,
                rho = if (.N >= 3 && uniqueN(gdi_centered) > 1 && uniqueN(ghi_centered) > 1)
                        suppressWarnings(cor(gdi_centered, ghi_centered, method = "spearman"))
                      else NA_real_), by = study][order(rho)]
for (i in seq_len(nrow(per)))
  log_msg("  %-24s n=%2d rho=%s", per$study[i], per$n[i],
          if (is.na(per$rho[i])) "  n/a" else sprintf("%+.3f", per$rho[i]))
okp <- per[!is.na(rho)]
if (nrow(okp)) {
  log_msg("")
  log_msg("  median %+.3f; %d of %d studies negative (%.0f%%)",
          stats::median(okp$rho), sum(okp$rho < 0), nrow(okp), 100*mean(okp$rho < 0))
  log_msg("  With 2-7 arms per study many of these are exactly +/-1. Read them")
  log_msg("  for pattern; the permutation test above is the inferential result.")
}

# =============================================================================
# STEP 5 - AUC AGAINST THE CURATED LABEL
# =============================================================================
log_msg("[5/7] AUC vs curated Challenge_Type")
auc_of <- function(g, pos, neg) {
  s1 <- g[pos & !is.na(g)]; s0 <- g[neg & !is.na(g)]
  if (length(s1) < 5 || length(s0) < 5) return(NA_real_)
  r <- rank(c(s1, s0))
  (sum(r[seq_along(s1)]) - length(s1)*(length(s1)+1)/2) / (length(s1)*length(s0))
}
ctrl_v <- pred$cls %in% HEALTHY; chal_v <- pred$cls %in% CHALLENGE
add_v  <- pred$cls %in% ADDITIVE
aucs <- list()
for (sp in c("centered", "absolute")) {
  g <- pred[[paste0("gdi_", sp)]]
  a_ch <- auc_of(g, chal_v, ctrl_v)  # challenged should score HIGHER (dysbiotic)
  a_ad <- auc_of(g, ctrl_v, add_v)   # control should score higher than additive
  aucs[[sp]] <- c(challenge = a_ch, additive = a_ad)
  log_msg("  %-9s challenged vs control : %s  (n=%d vs %d)", sp,
          if (is.na(a_ch)) "n/a" else sprintf("%.3f", a_ch), sum(chal_v), sum(ctrl_v))
  log_msg("  %-9s control vs additive   : %s  (n=%d vs %d)", sp,
          if (is.na(a_ad)) "n/a" else sprintf("%.3f", a_ad), sum(ctrl_v), sum(add_v))
}
log_msg("  0.5 = chance; >0.70 acceptable, >0.80 excellent.")
# This previously read that 04e validated the GHI for the additive contrast
# and found it null for challenge, so the GDI should be expected to inherit
# that asymmetry. That was true of an early run and has not been for some
# time: on the current seven-domain index every contrast in 04e is
# significant, including challenge (Combined +0.846, p = 0.0092, k = 24).
# The note is corrected rather than deleted, because the reasoning it
# contained still matters - a GDI trained on the GHI can only be as good as
# the target, and where the GDI fails the GHI's own performance is the first
# thing to check.
log_msg("  The GDI is trained on the GHI, so its ceiling is the target's own")
log_msg("  discrimination. On the current index that is AUC 0.651 across all")
log_msg("  scored samples and 0.884 for samples with three or more domains")
log_msg("  measured, so a weak GDI here is not explained by a weak target.")

# =============================================================================
# STEP 6 - PANEL CONSISTENCY
# =============================================================================
log_msg("[6/7] Panel consistency across folds")
tab <- sort(table(unlist(fold_panels)), decreasing = TRUE); nf <- length(fold_panels)
for (g in names(tab))
  log_msg("  %-40s %3d/%3d (%3.0f%%)%s", g, tab[[g]], nf, 100*tab[[g]]/nf,
          if (!is.null(m$panel) && g %in% m$panel) "  [in published panel]" else "")
core <- names(tab)[tab / nf >= 0.70]
log_msg("  appearing in >=70%% of folds: %d genera", length(core))
if (length(core)) log_msg("    %s", paste(core, collapse = ", "))
log_msg("  NOTE: fold appearance is not equivalent to the bootstrap stability")
log_msg("  reported in 06. Leave-one-out omits 1 of %d studies per iteration,", nf)
log_msg("  whereas stability selection omits 50%%, so fold appearance is")
log_msg("  systematically higher. The two are reported separately.")

# =============================================================================
# STEP 7 - WRITE
# =============================================================================
log_msg("[7/7] Writing")
gd <- file.path(R_OUT_DIR, "gdi")
fwrite(pred, file.path(gd, "loso_predictions.tsv"), sep = "\t")
fwrite(per,  file.path(gd, "loso_per_study.tsv"),   sep = "\t")

lines <- c("GDI LEAVE-ONE-STUDY-OUT VALIDATION",
  strrep("=", 74), "",
  "METHOD", strrep("-", 74),
  "Taxon selection and centroid estimation were performed within each training",
  "fold. The held-out study contributed to neither feature selection, centroid",
  "estimation, nor the median split defining the two reference groups.", "",
  sprintf("Arms scored                            : %d", nrow(pred)),
  sprintf("Studies held out                       : %d", uniqueN(pred$study)),
  sprintf("Folds using the reduced fallback panel : %d", n_fallback), "",
  "CORRELATION WITH THE GUT HEALTH INDEX", strrep("-", 74),
  "A negative coefficient is the expected direction, the GDI increasing with",
  "dysbiosis and the GHI with health. P-values are from permutation of the",
  "outcome within study, which preserves the clustering of arms within studies.",
  "")
for (sp in names(pool))
  lines <- c(lines, sprintf("  %-9s vs %-13s rho = %+.3f   p = %.4f",
                            sp, OUTCOME[[sp]], pool[[sp]]$rho, pool[[sp]]$p))
lines <- c(lines, "",
  "DISCRIMINATION BY INTERVENTION CLASS", strrep("-", 74),
  "Area under the ROC curve.", "")
for (sp in names(aucs))
  lines <- c(lines, sprintf("  %-9s challenged/control %s   control/additive %s", sp,
    ifelse(is.na(aucs[[sp]]["challenge"]), "n/a", sprintf("%.3f", aucs[[sp]]["challenge"])),
    ifelse(is.na(aucs[[sp]]["additive"]),  "n/a", sprintf("%.3f", aucs[[sp]]["additive"]))))
lines <- c(lines, "", "NOTES", strrep("-", 74),
  "1. Two scoring spaces are reported. In the centred space each arm is",
  "   expressed relative to the mean of its own study, which requires several",
  "   samples from one trial. In the absolute space arms are scored on raw CLR",
  "   values without a reference group. Each is evaluated against the outcome",
  "   on the corresponding scale.",
  "",
  "2. Taxon selection was performed on within-study centred data in both cases,",
  "   as the uncentred fit did not improve on the null model (cross-validated",
  "   error 151.96 against a null variance of 151.49).",
  "",
  "3. The differential abundance pre-filter was applied once to the full corpus",
  "   rather than within each fold. The candidate pool entering each fold was",
  "   therefore informed by the held-out study, and reported performance is to",
  "   that extent optimistic. All subsequent steps were fold-specific. The",
  "   pre-filter reduced 84 genera to 29, constraining rather than determining",
  "   the selected panel.",
  "",
  "4. Studies contribute between 2 and 7 arms. Per-study correlations are",
  "   reported for completeness but are unstable at this size and frequently",
  "   take the values +1 or -1. Inference rests on the pooled estimate.")
writeLines(lines, file.path(gd, "loso_summary.txt"))
log_msg("predictions -> %s", file.path(gd, "loso_predictions.tsv"))
log_msg("summary     -> %s", file.path(gd, "loso_summary.txt"))
log_msg("=== DONE ===")
