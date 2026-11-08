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
# Registered primary is "1se"; "min" is reported alongside as a sensitivity.
# Settable from the environment so the two can be compared without editing:
#   LAMBDA_RULE=min Rscript 06_lasso_selection.R
LAMBDA_RULE     <- { e <- Sys.getenv("LAMBDA_RULE", ""); if (nzchar(e)) e else "1se" }
if (!LAMBDA_RULE %in% c("1se","min")) stop("LAMBDA_RULE must be '1se' or 'min'")
set.seed(42)              # reproducibility

log_msg("[1/6] Loading MaAsLin2 output")
inp <- file.path(R_OUT_DIR, "gdi", "gdi_input.rds")
if (!file.exists(inp)) stop(sprintf("Not found: %s\nRun 05_maaslin2.R first.", inp))
d <- readRDS(inp)
clr <- d$clr; meta <- d$meta; selected <- d$selected
log_msg("  %d genera x %d runs; %d flagged by MaAsLin2", nrow(clr), ncol(clr), length(selected))
if (length(selected) < 2) stop("Fewer than 2 genera passed MaAsLin2 - nothing to select from.")

# ---- restrict to the MaAsLin2 survivors --------------------------------------
# BUGFIX (found in systematic audit, 13 Sep 2026): matching selected names
# against clr rownames by exact string equality silently drops any genus
# whose name MaAsLin2 sanitised differently from the CLR table. MaAsLin2
# replaces characters outside [A-Za-z0-9_] with periods when it writes
# feature names (e.g. "Escherichia-Shigella" -> "Escherichia.Shigella",
# "[Eubacterium]_coprostanoligenes_group" -> ".Eubacterium._coprostanoligenes_group"),
# while the CLR table keeps the original taxonomy-derived names with dashes
# and brackets intact. Of 28 genera MaAsLin2 flagged, exact-string matching
# found only 22 - six real, significant genera were silently excluded from
# every downstream LASSO, elastic net, random forest, percentile
# normalisation, LOSO validation, and axis analysis in this project, with no
# warning anywhere. Matching is now done on a canonicalised form (letters
# and digits only, case-insensitive) so sanitisation differences cannot
# cause a silent drop, and every match is logged explicitly.
canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
match_idx <- match(canon(selected), canon(rownames(clr)))
n_unmatched <- sum(is.na(match_idx))
if (n_unmatched > 0) {
  log_msg("  WARNING: %d MaAsLin2-flagged genus name(s) still have no match after",
          n_unmatched)
  log_msg("  canonicalisation - printing them so this cannot fail silently:")
  for (nm in selected[is.na(match_idx)]) log_msg("    unmatched: %s", nm)
}
keep <- rownames(clr) %in% rownames(clr)[na.omit(match_idx)]
clr <- clr[keep, , drop = FALSE]
log_msg("  restricted to %d genera (of %d flagged; %d matched by canonical name)",
        nrow(clr), length(selected), sum(!is.na(match_idx)))

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
# WITHIN-STUDY CENTERING
#
# Diagnosed directly, not assumed: fitting on raw arm-level X and y gave a
# cross-validated error WORSE than predicting the mean (151.96 vs 151.49 null
# variance), for every choice of fold structure - the model has nothing to
# select on the between-study scale. MaAsLin2, by contrast, found real
# associations for 28 genera, because it modelled study as a random effect,
# which is a WITHIN-study test: does this genus track GHI among the arms of
# the same study, after removing that study's own baseline level.
#
# Collapsing runs to arm means and fitting LASSO directly on those means
# discards exactly the within-study contrast MaAsLin2 used, and leaves mostly
# between-study variation - which reflects lab, breed, diet base and age far
# more than gut health. That is why the direct fit found nothing: it was
# being asked a different, harder question than the one with a real answer.
#
# The fix is to subtract each study's own mean from both X and y before
# fitting, so the model is forced to use within-study contrasts - the same
# information MaAsLin2 used, expressed as a single regression rather than
# per-genus mixed models. A study contributing only one arm carries no
# within-study contrast and is dropped, since centering a single point yields
# exactly zero with no information.
# -----------------------------------------------------------------------------
study_n <- table(study)
single_arm_studies <- names(study_n)[study_n < 2]
if (length(single_arm_studies) > 0) {
  log_msg("  %d single-arm stud%s dropped (no within-study contrast possible): %s",
          length(single_arm_studies), if (length(single_arm_studies)==1) "y" else "ies",
          paste(single_arm_studies, collapse = ", "))
  keep2 <- !(study %in% single_arm_studies)
  X <- X[keep2, , drop = FALSE]; y <- y[keep2]; study <- study[keep2]
}
log_msg("  %d arms across %d multi-arm studies remain for within-study fitting",
        nrow(X), length(unique(study)))

Xc <- X; yc <- y
for (s in unique(study)) {
  idx <- study == s
  Xc[idx, ] <- sweep(X[idx, , drop = FALSE], 2, colMeans(X[idx, , drop = FALSE]), "-")
  yc[idx]   <- y[idx] - mean(y[idx])
}
log_msg("  centered: y SD before=%.3f, within-study SD after=%.3f", stats::sd(y), stats::sd(yc))

n <- nrow(Xc)
n_folds <- max(3, min(10, floor(n / 8)))
ord <- order(yc)
fold_id <- integer(n)
fold_id[ord] <- ((seq_len(n) - 1) %% n_folds) + 1

cvf <- cv.glmnet(Xc, yc, alpha = 1, foldid = fold_id, standardize = TRUE,
                 intercept = FALSE)   # centered data has no intercept to fit
lam <- if (LAMBDA_RULE == "1se") cvf$lambda.1se else cvf$lambda.min

null_var <- stats::var(yc)
log_msg("  lambda.min CV error=%.4f vs within-study null variance=%.4f (ratio %.3f)",
        min(cvf$cvm), null_var, min(cvf$cvm) / null_var)
if (min(cvf$cvm) / null_var < 0.98) {
  log_msg("  the model beats the within-study mean - there is signal to select from")
} else {
  log_msg("  *** the model still does not beat the within-study mean.")
  log_msg("  *** This would mean the MaAsLin2 associations do not survive an")
  log_msg("  *** independent multi-genus fit, even within-study - a genuinely")
  log_msg("  *** informative negative result, not a coding problem.")
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
counts <- setNames(numeric(ncol(Xc)), colnames(Xc))
n_ok <- 0
for (b in seq_len(N_BOOT)) {
  take <- sample(studies, n_take)
  idx <- which(study %in% take)
  if (length(idx) < 10 || length(unique(yc[idx])) < 5) next
  fit <- tryCatch(
    cv.glmnet(Xc[idx, , drop = FALSE], yc[idx], alpha = 1,
              nfolds = min(5, length(idx)), standardize = TRUE, intercept = FALSE),
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
pred <- as.numeric(predict(cvf, newx = Xc, s = lam))
r2 <- 1 - sum((yc - pred)^2) / sum((yc - mean(yc))^2)
log_msg("  (within-study centered scale)")
log_msg("  R2=%.3f  Pearson r=%.3f  Spearman rho=%.3f",
        r2, stats::cor(yc, pred), stats::cor(yc, pred, method = "spearman"))
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
# coef_1se is kept under its own name for continuity. coef_sel holds the
# coefficients at whichever lambda this run used, and is what 07 reads, so the
# panel rule follows LAMBDA_RULE rather than always following lambda.1se.
sel_used <- if (LAMBDA_RULE == "1se") sel_1se else sel_min
fwrite(data.table(genus = names(freq), stability = as.numeric(freq),
                  coef_1se = as.numeric(sapply(names(freq), function(g)
                    if (g %in% names(sel_1se)) sel_1se[[g]] else 0)),
                  coef_sel = as.numeric(sapply(names(freq), function(g)
                    if (g %in% names(sel_used)) sel_used[[g]] else 0)),
                  lambda_rule = LAMBDA_RULE
       )[order(-stability)], file.path(gd, "gdi_stability.tsv"), sep = "\t")
saveRDS(list(X = X, y = y, Xc = Xc, yc = yc, study = study, panel = panel, freq = freq,
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
