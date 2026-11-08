#!/usr/bin/env Rscript
# =============================================================================
# 03_batch_correction.R
#
# Corrects for study-to-study (batch) technical variation in the genus-level
# feature table, then applies a CLR transformation to the corrected output.
#
# HOW TO RUN: from a terminal (Rscript 03_batch_correction.R) or inside
# RStudio (open the file, click Source). No path editing required; see
# 00_r_config.R for the expected repository layout. Must be run AFTER
# 01_genus_agglomeration.R, since it consumes that script's output.
#
# STEP ORDER (verified against MMUPHin's own documentation, source code, and
# multiple published applications; this order is a hard requirement, not a
# style choice):
#
#   1. Rare-taxon filtering        (already done in 01_genus_agglomeration.R)
#   2. Convert counts -> proportions
#   3. adjust_batch (MMUPHin)      runs on proportions/counts
#   4. Zero-replacement (CZM)      on the BATCH-CORRECTED proportions
#   5. CLR transformation          on the zero-replaced, batch-corrected data
#   6. Export                      CLR matrix used by all downstream steps
#
# WHY THIS ORDER, SPECIFICALLY: MMUPHin's adjust_batch performs a zero-
# inflated empirical Bayes correction designed for non-negative compositional
# data (proportions or counts), and distinguishes true structural zeros from
# sampling zeros as part of that model. It returns output "normalized to the
# same per-sample total abundance as feature_abd" - a statement that is only
# meaningful for compositional data, not CLR values (which are already
# log-ratios, centered at zero, and can be negative). Passing CLR-transformed
# data into adjust_batch is reported to fail outright ("Feature table must be
# non-negative for normalization"), and no published MMUPHin-based
# meta-analysis surveyed applies CLR before batch correction; CLR is
# consistently applied afterward, to the corrected output, before PCA or
# differential abundance testing.
#
# COVARIATE PROTECTION: adjust_batch is told which biological variable(s) to
# preserve while removing the batch (Study_ID) effect, via the `covariates`
# argument. This choice is informed directly by 02b_batch_confounding_check.R:
#   - Intervention_Class is NOT strongly confounded with Study_ID (in the
#     reference corpus, the large majority of studies contained more than one
#     Intervention_Class level) and is safe to protect.
#   - Segment IS strongly confounded with Study_ID (in the reference corpus,
#     the large majority of studies sampled only one gut segment). Passing
#     Segment as a protected covariate would ask adjust_batch to separate two
#     variables that are nearly collinear in this data, which the model
#     cannot reliably do. Segment is deliberately NOT protected here; any
#     segment-stratified comparison downstream should be treated as
#     confounded with batch and interpreted with that caveat, not corrected
#     for as though it were independent of study identity.
# Re-run 02b_batch_confounding_check.R against the current metadata before
# changing PROTECTED_COVARIATES, and update this comment to match.
#
# WHAT HAS AND HAS NOT BEEN TESTED: the data loading, sample-alignment, and
# CLR math in this script have been verified against synthetic data. The
# MMUPHin::adjust_batch() call itself requires Bioconductor packages
# (MMUPHin, zCompositions) not available in every environment; it must be
# verified on first real run against actual data, per the checklist printed
# at the end of this script.
# =============================================================================

# ---- Locate this script's folder and load shared configuration -------------
.locate_and_source_config <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    this_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else if (requireNamespace("rstudioapi", quietly = TRUE) &&
             isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    this_dir <- if (!is.null(ctx) && nzchar(ctx$path)) dirname(normalizePath(ctx$path)) else normalizePath(".")
  } else {
    this_dir <- normalizePath(".")
  }
  source(file.path(this_dir, "00_r_config.R"))
}
.locate_and_source_config()

# ---- Additional packages beyond 00_r_config.R's base check -----------------
# MMUPHin and zCompositions are Bioconductor packages, installed via
# BiocManager rather than install.packages(), so they are checked separately
# from the CRAN packages in 00_r_config.R.
bioc_packages <- c("MMUPHin", "zCompositions")
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  stop(sprintf(
    paste0("Missing required Bioconductor packages: %s\n",
           "Install with:\n",
           "  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")\n",
           "  BiocManager::install(c(%s))"),
    paste(missing_bioc, collapse = ", "),
    paste(sprintf('"%s"', missing_bioc), collapse = ", ")
  ))
}
suppressPackageStartupMessages({
  library(MMUPHin)
  library(zCompositions)
})

# ---- Parameters --------------------------------------------------------------
BATCH_VAR <- "Study_ID"
# See the header comment above for why Intervention_Class is protected and
# Segment is not. Re-derive this from 02b's output before changing it.
#
# AGE: Age_d is NOT a column in the QIIME2 sample metadata (it lives only in
# the biomarker workbook, joined in below by Study_ID + Treatment_Group). The
# registered pipeline's own variance-partition test found age explains far
# more community variance than anything else measured (R2 0.12), so leaving
# it unprotected risks adjust_batch treating real age-driven signal as batch
# noise wherever age happens to correlate with study - which it often does,
# since many studies sample at a single age. Protecting it tells adjust_batch
# to preserve age-associated variation while still removing the Study_ID
# batch effect.
#
# AGE_COVARIATE=no restores the original single-covariate correction, so the
# two can be compared as a sensitivity analysis.
AGE_COVARIATE <- { e <- tolower(Sys.getenv("AGE_COVARIATE", "")); !(e %in% c("no","false","0")) }
PROTECTED_COVARIATES <- c("Batch_Covariate")   # Age_d appended below if used and available

# -----------------------------------------------------------------------------
# STEP 1 - LOAD THE GENUS-LEVEL COUNT TABLE (output of 01_genus_agglomeration.R)
# -----------------------------------------------------------------------------
genus_counts_path <- file.path(R_OUT_DIR, "genus", "genus_counts.tsv")
if (!file.exists(genus_counts_path)) {
  stop(sprintf("Genus count table not found at %s.\nRun 01_genus_agglomeration.R first.",
               genus_counts_path))
}
log_msg("[1/6] Loading genus count table from %s", genus_counts_path)
genus_dt <- fread(genus_counts_path, sep = "\t", header = TRUE)

genus_names <- genus_dt$genus
count_mat <- as.matrix(genus_dt[, -"genus"])
rownames(count_mat) <- genus_names
# count_mat is now features (genera) x samples, matching MMUPHin's expected
# feature-by-sample orientation directly - no transposition needed.

log_msg("  %d genera x %d samples", nrow(count_mat), ncol(count_mat))

# -----------------------------------------------------------------------------
# STEP 2 - LOAD METADATA AND ALIGN SAMPLES
# adjust_batch requires the feature matrix's column names and the metadata
# data frame's row names to refer to the same samples, in a data frame
# (not data.table) with rownames set explicitly. Samples present in one
# input but not the other are dropped, and the exact count dropped is
# reported - a silent mismatch here would misattribute batch/covariate
# labels to the wrong samples.
# -----------------------------------------------------------------------------
log_msg("[2/6] Loading and aligning metadata from %s", METADATA_TSV)
meta_raw <- fread(METADATA_TSV, sep = "\t", header = FALSE, skip = 2, na.strings = c("", "NA"))
meta_header <- fread(METADATA_TSV, sep = "\t", header = FALSE, nrows = 1) %>% unlist(use.names = FALSE)
setnames(meta_raw, meta_header)

required_meta_cols <- c("sample-id", BATCH_VAR, PROTECTED_COVARIATES)
missing_cols <- setdiff(required_meta_cols, names(meta_raw))
if (length(missing_cols) > 0) {
  stop(sprintf("Metadata is missing required column(s): %s", paste(missing_cols, collapse = ", ")))
}

meta_df <- as.data.frame(meta_raw[, ..required_meta_cols])
rownames(meta_df) <- meta_df[["sample-id"]]
meta_df[["sample-id"]] <- NULL

common_samples <- intersect(colnames(count_mat), rownames(meta_df))
n_dropped_from_counts <- ncol(count_mat) - length(common_samples)
n_dropped_from_meta   <- nrow(meta_df) - length(common_samples)
if (n_dropped_from_counts > 0) {
  log_msg("  WARNING: %d samples in the feature table have no metadata row and were dropped", n_dropped_from_counts)
}
if (n_dropped_from_meta > 0) {
  log_msg("  NOTE: %d metadata rows have no matching feature-table sample (expected - not all samples survive QIIME 2 filtering)", n_dropped_from_meta)
}

count_mat <- count_mat[, common_samples, drop = FALSE]
meta_df   <- meta_df[common_samples, , drop = FALSE]
meta_df[[BATCH_VAR]] <- as.factor(meta_df[[BATCH_VAR]])

log_msg("  %d samples retained after alignment", length(common_samples))
log_msg("  %d batches (studies) present", nlevels(meta_df[[BATCH_VAR]]))

# -----------------------------------------------------------------------------
# JOIN AGE FROM THE BIOMARKER WORKBOOK
#
# Age_d does not exist in the QIIME2 sample metadata, so it is joined in here
# by Study_ID + Treatment_Group, the same key used throughout the pipeline to
# connect sequence data to biomarker data. A run whose arm has more than one
# recorded age (a study sampling at multiple timepoints under one label) gets
# the mean of those ages - the same handling as any other duplicate-key
# collapse elsewhere in this project.
#
# adjust_batch cannot accept a covariate with missing values, so a run missing
# age is EXCLUDED from the batch correction entirely when AGE_COVARIATE is
# active, rather than imputed. Imputing a covariate risks manufacturing the
# very structure adjust_batch is being asked to preserve; dropping and
# reporting the count is the more conservative choice.
# -----------------------------------------------------------------------------
if (AGE_COVARIATE) {
  if (!file.exists(METADATA_XLSX)) {
    log_msg("  AGE_COVARIATE requested but %s not found - proceeding without age", METADATA_XLSX)
  } else {
    bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                         stringsAsFactors = FALSE)
    if (!"Age_d" %in% names(bio)) {
      log_msg("  AGE_COVARIATE requested but Age_d not found in %s - proceeding without age",
              BIOMARKER_SHEET)
    } else {
      bk <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||")
      age_lut <- tapply(suppressWarnings(as.numeric(bio$Age_d)), bk, mean, na.rm = TRUE)
      rk <- paste(trimws(meta_df[[BATCH_VAR]]), trimws(meta_df$Intervention_Class), sep = "||")
      # Intervention_Class is not the same key as Treatment_Group; fall back to
      # a direct Treatment_Group column if present in meta_df, else skip.
      if ("Treatment_Group" %in% names(meta_raw)) {
        rk <- paste(trimws(meta_raw$Study_ID[match(rownames(meta_df), meta_raw[["sample-id"]])]),
                   trimws(meta_raw$Treatment_Group[match(rownames(meta_df), meta_raw[["sample-id"]])]),
                   sep = "||")
      }
      meta_df$Age_d <- unname(age_lut[rk])
      n_age_missing <- sum(is.na(meta_df$Age_d))
      log_msg("  age joined: %d of %d samples have a usable age (%d missing, will be excluded)",
              sum(!is.na(meta_df$Age_d)), nrow(meta_df), n_age_missing)
      if (sum(!is.na(meta_df$Age_d)) >= 30) {
        keep_age <- !is.na(meta_df$Age_d)
        count_mat <- count_mat[, keep_age, drop = FALSE]
        meta_df   <- meta_df[keep_age, , drop = FALSE]
        PROTECTED_COVARIATES <- c(PROTECTED_COVARIATES, "Age_d")
        log_msg("  Age_d added to protected covariates; %d samples retained", nrow(meta_df))
      } else {
        log_msg("  fewer than 30 samples have a usable age - proceeding WITHOUT the age covariate")
      }
    }
  }
}

# -----------------------------------------------------------------------------
# STEP 3 - COUNTS -> PROPORTIONS
# adjust_batch accepts proportions or counts per its documentation, but the
# package's own worked example uses proportions specifically, so this script
# matches that more thoroughly tested path rather than passing raw counts.
# -----------------------------------------------------------------------------
log_msg("[3/6] Converting counts to relative abundance (proportions)")
prop_mat <- sweep(count_mat, 2, colSums(count_mat), FUN = "/")

# A sample with zero total reads after filtering would produce NaN here and
# must be caught explicitly rather than silently propagating into adjust_batch.
zero_total_samples <- colnames(count_mat)[colSums(count_mat) == 0]
if (length(zero_total_samples) > 0) {
  stop(sprintf(
    "%d sample(s) have zero total genus-level counts after filtering: %s\nThese must be removed before batch correction.",
    length(zero_total_samples), paste(head(zero_total_samples, 5), collapse = ", ")
  ))
}

# -----------------------------------------------------------------------------
# STEP 4 - MMUPHin BATCH CORRECTION
# Runs on proportions, protecting only the covariate(s) confirmed not to be
# collinear with batch (see header comment). Output is a proportions matrix
# on the same scale as the input.
# -----------------------------------------------------------------------------
log_msg("[4/6] Running MMUPHin adjust_batch (batch = %s, protected covariate(s) = %s)",
        BATCH_VAR, paste(PROTECTED_COVARIATES, collapse = ", "))

fit_adjust_batch <- adjust_batch(
  feature_abd = prop_mat,
  batch       = BATCH_VAR,
  covariates  = PROTECTED_COVARIATES,
  data        = meta_df
)
prop_mat_adj <- fit_adjust_batch$feature_abd_adj

out_diag <- file.path(R_OUT_DIR, "genus", "mmuphin_diagnostic.pdf")
log_msg("  batch-corrected proportions produced (%d x %d)", nrow(prop_mat_adj), ncol(prop_mat_adj))

# -----------------------------------------------------------------------------
# STEP 5 - ZERO-REPLACEMENT (CZM) THEN CLR, ON THE BATCH-CORRECTED OUTPUT
# CLR requires strictly positive values (it is a log-ratio), but compositional
# proportions legitimately contain exact zeros (a genus absent from a
# sample). cmultRepl's count zero multiplicative (CZM) method replaces zeros
# with small positive values while preserving the compositional closure
# (each sample's proportions still sum to 1) - this is standard practice in
# compositional data analysis and is applied here, not before batch
# correction, so that the zero model MMUPHin itself uses is not disturbed.
#
# IMPORTANT - cmultRepl SILENTLY DROPS rows/columns that are >80% zeros
# (controlled by its z.delete argument, default TRUE), because zero
# replacement is not statistically reliable with that little non-zero data
# to estimate from. Left unchecked this can shrink the feature/sample set
# without any record of what was removed or why. This step explicitly
# captures the before/after names on both axes and reports the exact
# identities dropped - never just a count - so a silent data loss cannot
# pass through unnoticed.
# -----------------------------------------------------------------------------
log_msg("[5/6] Zero-replacement (CZM) and CLR transformation")

# cmultRepl expects samples as ROWS, features as COLUMNS - the transpose of
# the feature-by-sample convention used everywhere else in this script.
prop_mat_adj_t <- t(prop_mat_adj)
genera_before   <- colnames(prop_mat_adj_t)
samples_before  <- rownames(prop_mat_adj_t)

prop_repl_t <- cmultRepl(prop_mat_adj_t, method = "CZM", output = "p-counts")

genera_after  <- colnames(prop_repl_t)
samples_after <- rownames(prop_repl_t)
genera_dropped_by_repl  <- setdiff(genera_before, genera_after)
samples_dropped_by_repl <- setdiff(samples_before, samples_after)

# Reported unconditionally (even when the dropped set is empty) so the
# summary file always states explicitly what happened here, rather than
# only mentioning it when something goes wrong.
log_msg("  cmultRepl: %d of %d genera dropped (>80%% zeros, cannot be zero-replaced reliably)",
        length(genera_dropped_by_repl), length(genera_before))
log_msg("  cmultRepl: %d of %d samples dropped (>80%% zeros, cannot be zero-replaced reliably)",
        length(samples_dropped_by_repl), length(samples_before))

if (length(genera_dropped_by_repl) > 0) {
  out_dropped_genera <- file.path(R_OUT_DIR, "genus", "cmultRepl_dropped_genera.txt")
  writeLines(genera_dropped_by_repl, out_dropped_genera)
  log_msg("  dropped genus names -> %s", out_dropped_genera)
}
if (length(samples_dropped_by_repl) > 0) {
  out_dropped_samples <- file.path(R_OUT_DIR, "genus", "cmultRepl_dropped_samples.txt")
  writeLines(samples_dropped_by_repl, out_dropped_samples)
  log_msg("  dropped sample IDs -> %s", out_dropped_samples)
}

# A large fraction dropped indicates the prevalence filter in
# 01_genus_agglomeration.R (currently 10%) is not strict enough relative to
# what cmultRepl itself requires (roughly >=20% non-zero, i.e. <80% zero) -
# the two filters are not currently aligned, and this warns explicitly
# rather than letting a large, unreviewed loss pass silently.
pct_genera_dropped <- 100 * length(genera_dropped_by_repl) / length(genera_before)
if (pct_genera_dropped > 10) {
  log_msg("  *** WARNING: %.1f%% of genera were dropped at the zero-replacement step.", pct_genera_dropped)
  log_msg("  *** This suggests the prevalence filter in 01_genus_agglomeration.R is looser")
  log_msg("  *** than what cmultRepl's 80%% zero-deletion threshold requires (roughly")
  log_msg("  *** equivalent to a >=20%% prevalence filter). Consider raising the prevalence")
  log_msg("  *** threshold there to reduce this loss, or confirm with your statistician")
  log_msg("  *** that this loss is acceptable before proceeding.")
}

prop_repl_t <- prop_repl_t / rowSums(prop_repl_t)   # renormalise to proportions after replacement

clr_transform <- function(row) {
  log(row) - mean(log(row))
}
clr_mat_t <- t(apply(prop_repl_t, 1, clr_transform))
clr_mat <- t(clr_mat_t)   # back to features x samples for consistency with count_mat/prop_mat

log_msg("  CLR matrix: %d genera x %d samples", nrow(clr_mat), ncol(clr_mat))

# -----------------------------------------------------------------------------
# STEP 6 - EXPORT
# Both the batch-corrected proportions and the final CLR matrix are written.
# The CLR matrix is what downstream PCA/GHI/MaAsLin2 steps should consume;
# the proportions matrix is retained for any step that specifically needs
# compositional (non-CLR) data, or for diagnostic comparison to the
# pre-correction proportions.
# -----------------------------------------------------------------------------
log_msg("[6/6] Writing output")

out_prop_adj <- file.path(R_OUT_DIR, "genus", "genus_proportions_batch_adjusted.tsv")
out_clr      <- file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")

fwrite(data.table(genus = rownames(prop_mat_adj), prop_mat_adj), out_prop_adj, sep = "\t")
fwrite(data.table(genus = rownames(clr_mat), clr_mat), out_clr, sep = "\t")

log_msg("  batch-adjusted proportions -> %s", out_prop_adj)
log_msg("  batch-adjusted CLR matrix  -> %s", out_clr)

summary_lines <- c(
  "=== Batch correction summary ===",
  sprintf("Samples: %d", ncol(count_mat)),
  sprintf("Genera: %d", nrow(count_mat)),
  sprintf("Batches (studies): %d", nlevels(meta_df[[BATCH_VAR]])),
  sprintf("Protected covariate(s): %s", paste(PROTECTED_COVARIATES, collapse = ", ")),
  "Segment was NOT protected (confounded with Study_ID; see 02b_batch_confounding_check.R).",
  "",
  "=== First-run verification checklist (check these manually) ===",
  "[ ] adjust_batch ran without error and without an unexpected warning about",
  "    the covariate model rank or convergence.",
  "[ ] fit_adjust_batch$batch_diagnostic (or the diagnostic PDF, if generated)",
  "    shows reduced study clustering relative to the pre-correction PCoA from",
  "    03a_merge_and_pcoa.sh.",
  "[ ] Re-running the Study_ID PERMANOVA on the batch-adjusted data (not yet",
  "    scripted here) shows a materially lower pseudo-F than the pre-correction",
  "    value, while Intervention_Class retains a detectable signal.",
  "[ ] No sample or genus was silently dropped beyond what this log reports."
)
out_summary <- file.path(R_OUT_DIR, "genus", "batch_correction_summary.txt")
writeLines(summary_lines, out_summary)
log_msg("  summary and verification checklist -> %s", out_summary)

log_msg("=== DONE ===")
