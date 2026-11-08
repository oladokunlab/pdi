#!/usr/bin/env Rscript
# =============================================================================
# 04f_convergent_validity.R
#
# Do the Morphology and Biomarker arms of the GHI measure the same thing?
#
# RUN:  Rscript 04f_convergent_validity.R
# TIME: seconds. Reads ghi_scores_pca.tsv only. Changes nothing.
#
# -----------------------------------------------------------------------------
# WHY THIS EXISTS, AND WHAT IT IS NOT
# -----------------------------------------------------------------------------
# GHI_Combined is the pre-specified primary index and is built as an arithmetic
# mean across domain scores. That mean is only a meaningful summary if the
# domains track a shared underlying construct. If the morphological arm and the
# biomarker arm are uncorrelated on samples that carry both, then "Combined" is
# not a composite of gut health - it is the average of two unrelated
# measurements, and every downstream model built against it inherits that.
#
# This is CONVERGENT VALIDITY and it is SECONDARY evidence. It is not a
# validation of the index:
#
#   - Morphology is not a gold standard. It is another imperfect proxy. If the
#     two arms disagree, this test cannot say which one is wrong.
#   - It is partly circular for GHI_Combined, which contains both arms.
#   - Primary validation remains the mixed model in 04e, which tests against
#     an EXTERNAL label (known intervention direction) rather than against
#     another index.
#
# The question here is deliberately coarse: correlated, or not correlated?
# The overlap is small, so the point estimate is imprecise and the confidence
# interval is the thing to read. Do not over-interpret the second decimal.
#
# -----------------------------------------------------------------------------
# HOW TO READ THE RESULT
# -----------------------------------------------------------------------------
#   CI excludes 0, r positive     -> the arms converge. Combined is coherent.
#                                    Proceed to the immunoglobulin split.
#   CI includes 0, or r near 0    -> STOP. Combined is averaging unrelated
#                                    signals. Reconsider it as the GDI target
#                                    before any LASSO work.
#   r strongly negative           -> a sign/orientation error somewhere upstream.
#                                    Check the direction-of-effect table in 04b.
#
# Both raw and 0-100 scores are reported. Pearson r is invariant to the linear
# rescale, so the two must agree; if they do not, something is wrong with the
# rescale and that is worth knowing.
#
# Within-study correlation is also reported. The pooled figure can be inflated
# purely by between-study differences - the same artefact that produced the
# R2 = 0.000 bug in the LASSO. The within-study number is the honest one.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

SCRIPT_DIR <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())

source(file.path(SCRIPT_DIR, "00_r_config.R"))

GHI_TSV <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
if (!file.exists(GHI_TSV))
  stop("Not found: ", GHI_TSV, "\n  Run 04b_ghi_pca.R first.")

ghi <- fread(GHI_TSV)
log_msg("[1/4] Loaded %d rows from %s", nrow(ghi), basename(GHI_TSV))

# --- overlap -----------------------------------------------------------------
# Samples carrying BOTH arms. Everything below is computed on these only.
both <- ghi[!is.na(GHI_Morphology_raw) & !is.na(GHI_Biomarker_raw)]

log_msg("[2/4] Overlap")
log_msg("  Morphology available : %d", sum(!is.na(ghi$GHI_Morphology_raw)))
log_msg("  Biomarker available  : %d", sum(!is.na(ghi$GHI_Biomarker_raw)))
log_msg("  BOTH (analysed here) : %d", nrow(both))
log_msg("  studies represented  : %d", uniqueN(both$Study_ID))

if (nrow(both) < 10)
  stop("Overlap too small (", nrow(both), ") for a meaningful correlation.")

# --- pooled correlation ------------------------------------------------------
# Pearson on raw scores is the headline number. Spearman is reported alongside
# because PC1-based scores are not guaranteed to relate linearly, and a large
# gap between the two is itself informative (outlier-driven vs monotone).
log_msg("[3/4] Correlation")

ct_raw  <- cor.test(both$GHI_Morphology_raw, both$GHI_Biomarker_raw,
                    method = "pearson")
ct_100  <- cor.test(both$GHI_Morphology,     both$GHI_Biomarker,
                    method = "pearson")
ct_sp   <- suppressWarnings(
             cor.test(both$GHI_Morphology_raw, both$GHI_Biomarker_raw,
                      method = "spearman"))

fmt_ct <- function(ct) {
  ci <- if (!is.null(ct$conf.int)) sprintf("[%+.3f, %+.3f]", ct$conf.int[1],
                                           ct$conf.int[2]) else "     -        "
  sprintf("r = %+.3f  95%% CI %s  p = %.4g", unname(ct$estimate), ci, ct$p.value)
}

log_msg("  raw scores    %s", fmt_ct(ct_raw))
log_msg("  0-100 scores  %s", fmt_ct(ct_100))
log_msg("  Spearman      %s", fmt_ct(ct_sp))

if (abs(unname(ct_raw$estimate) - unname(ct_100$estimate)) > 1e-6)
  log_msg("  WARNING: raw and rescaled r differ. The 0-100 rescale should be")
log_msg("           linear and cannot change Pearson r. Check step 6 of 04b.")

# --- within-study correlation ------------------------------------------------
# Pooling across studies lets between-study variation masquerade as agreement.
# Centring each arm within its own study removes that, leaving only the
# within-study contrast - the same correction that fixed the LASSO.
# Studies contributing a single overlapping sample carry no within-study
# contrast and drop out.
ws <- copy(both)
ws[, n_study := .N, by = Study_ID]
ws <- ws[n_study >= 3]

within_r <- NA_real_
if (nrow(ws) >= 10 && uniqueN(ws$Study_ID) >= 2) {
  ws[, `:=`(m_c = GHI_Morphology_raw - mean(GHI_Morphology_raw),
            b_c = GHI_Biomarker_raw  - mean(GHI_Biomarker_raw)), by = Study_ID]
  ct_w <- cor.test(ws$m_c, ws$b_c, method = "pearson")
  within_r <- unname(ct_w$estimate)
  log_msg("  within-study  %s", fmt_ct(ct_w))
  log_msg("                (%d samples, %d studies with >=3 overlapping)",
          nrow(ws), uniqueN(ws$Study_ID))
} else {
  log_msg("  within-study  not computable - too few studies with >=3 overlap")
}

# --- verdict -----------------------------------------------------------------
log_msg("[4/4] Verdict")

r  <- unname(ct_raw$estimate)
lo <- ct_raw$conf.int[1]
hi <- ct_raw$conf.int[2]

verdict <- if (lo > 0 && r >= 0.30) {
  "CONVERGE - Combined is coherent. Proceed."
} else if (lo > 0) {
  "WEAK BUT POSITIVE - Combined defensible, state the weak convergence."
} else if (hi < 0) {
  "NEGATIVE - suspect a direction-of-effect error. Check 04b DOMAINS table."
} else {
  "NOT ESTABLISHED - CI includes zero. Reconsider Combined as the GDI target."
}
log_msg("  %s", verdict)

if (!is.na(within_r) && lo > 0 && within_r < 0.15)
  log_msg("  CAUTION: pooled agreement is not reproduced within studies (%+.3f).",
          within_r)

# --- write -------------------------------------------------------------------
lines <- c(
  "GHI CONVERGENT VALIDITY - Morphology vs Biomarker",
  "Secondary evidence. Not a validation of the index; see 04e for primary.",
  strrep("=", 74), "",
  sprintf("Samples with both arms : %d", nrow(both)),
  sprintf("Studies represented    : %d", uniqueN(both$Study_ID)), "",
  sprintf("Pearson  (raw)   : %s", fmt_ct(ct_raw)),
  sprintf("Pearson  (0-100) : %s", fmt_ct(ct_100)),
  sprintf("Spearman (raw)   : %s", fmt_ct(ct_sp)),
  sprintf("Within-study     : %s",
          if (is.na(within_r)) "not computable" else sprintf("%+.3f", within_r)),
  "",
  sprintf("VERDICT: %s", verdict),
  "",
  "Morphology is not a gold standard. A disagreement here does not identify",
  "which arm is at fault; it only shows the two are not measuring one thing.",
  "The overlap is small, so read the interval rather than the point estimate."
)
outs <- file.path(R_OUT_DIR, "ghi", "ghi_convergent_validity.txt")
writeLines(lines, outs)
log_msg("summary -> %s", outs)

outd <- file.path(R_OUT_DIR, "ghi", "ghi_convergent_pairs.tsv")
fwrite(both[, intersect(c("Study_ID","Sample_ID","Treatment_Group",
                          "GHI_Morphology_raw","GHI_Biomarker_raw",
                          "GHI_Morphology","GHI_Biomarker","GHI_n_domains"),
                        names(both)), with = FALSE], outd, sep = "\t")
log_msg("paired scores -> %s", outd)
