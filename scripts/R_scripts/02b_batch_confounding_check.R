#!/usr/bin/env Rscript
# =============================================================================
# 02b_batch_confounding_check.R
#
# Checks whether Study_ID (the batch-correction variable) is confounded with
# biologically meaningful covariates, before batch correction is applied.
#
# HOW TO RUN: from a terminal (Rscript 02b_batch_confounding_check.R) or
# inside RStudio (open the file, click Source). No path editing required;
# see 00_r_config.R for the expected repository layout.
#
# WHY THIS RUNS BEFORE BATCH CORRECTION: batch correction removes variance
# associated with the batch (study) variable. If a covariate of interest -
# segment, intervention class, bird type - is itself unevenly distributed
# across studies, batch correction risks removing genuine biological signal
# along with the technical batch effect, because the two are statistically
# entangled and the correction cannot distinguish them.
#
# OUTPUT: for each categorical covariate, a contingency table against
# Study_ID plus a chi-squared test of independence; for each numerical
# covariate, per-study quartiles plus a Kruskal-Wallis test. A covariate
# strongly associated with Study_ID (visually uneven contingency table, or a
# significant result concentrated in a few studies) is flagged for review
# before batch correction proceeds.
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

# -----------------------------------------------------------------------------
# STEP 1 - LOAD METADATA
# -----------------------------------------------------------------------------
log_msg("[1/4] Loading metadata from %s", METADATA_TSV)

meta_raw <- fread(METADATA_TSV, sep = "\t", header = FALSE, skip = 2, na.strings = c("", "NA"))
meta_header <- fread(METADATA_TSV, sep = "\t", header = FALSE, nrows = 1) %>% unlist(use.names = FALSE)
setnames(meta_raw, meta_header)

log_msg("  %d samples, %d metadata columns", nrow(meta_raw), ncol(meta_raw))

# -----------------------------------------------------------------------------
# STEP 2 - IDENTIFY COVARIATES TO CHECK
# Categorical: variables describing intervention, segment, or design that a
# downstream model will use as a covariate or grouping variable.
# Numerical: variables where an uneven distribution across studies (e.g. one
# study measured entirely at one age) would indicate a design confound.
# Only columns actually present in the metadata file are checked, so this
# script runs unmodified even if some optional columns are absent.
# -----------------------------------------------------------------------------
categorical_covariates <- intersect(
  c("Segment", "Intervention_Class", "Challenge_Subtype", "Bird_Type",
    "Segment_matched", "PCoA_Group"),
  names(meta_raw)
)
numerical_covariates <- intersect(
  c("Read_count", "Birds_per_sample"),
  names(meta_raw)
)

log_msg("[2/4] Checking %d categorical and %d numerical covariates against Study_ID",
        length(categorical_covariates), length(numerical_covariates))

# -----------------------------------------------------------------------------
# STEP 3 - CATEGORICAL COVARIATES: CONTINGENCY TABLE + CHI-SQUARED
# withCallingHandlers is used rather than a plain tryCatch(..., warning=)
# handler, because the latter would discard the test RESULT entirely
# whenever chisq.test's own "approximation may be incorrect" warning fires -
# which is precisely the informative case this check needs to capture and
# report, not discard. A warning here is recorded as a flag on the result,
# not treated as a failure.
# -----------------------------------------------------------------------------
cat_results <- list()
for (cov in categorical_covariates) {
  tab <- table(meta_raw$Study_ID, meta_raw[[cov]], useNA = "ifany")
  out_path <- file.path(R_OUT_DIR, "confounding", sprintf("contingency_%s.tsv", cov))
  write.table(as.data.frame.matrix(tab), out_path, sep = "\t", quote = FALSE, col.names = NA)

  had_low_count_warning <- FALSE
  test <- tryCatch(
    withCallingHandlers(
      chisq.test(tab),
      warning = function(w) {
        if (grepl("Chi-squared approximation", conditionMessage(w))) {
          had_low_count_warning <<- TRUE
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) NULL
  )
  pct_low_expected <- if (!is.null(test)) mean(test$expected < 5) * 100 else NA

  cat_results[[cov]] <- list(
    covariate = cov,
    p_value = if (!is.null(test)) test$p.value else NA,
    pct_cells_low_expected = pct_low_expected
  )
  log_msg("  %-22s p=%s  (%.0f%% of cells have expected count <5%s)",
          cov,
          if (!is.null(test)) sprintf("%.4f", test$p.value) else "NA (test failed)",
          ifelse(is.na(pct_low_expected), 0, pct_low_expected),
          ifelse(!is.na(pct_low_expected) && pct_low_expected > 20,
                 " - chi-squared unreliable here, read the contingency table directly", ""))
}

# -----------------------------------------------------------------------------
# STEP 4 - NUMERICAL COVARIATES: PER-STUDY DISTRIBUTION
# Reports quartiles per study rather than only a single test statistic,
# since the object of interest is which studies deviate, not just whether
# any do.
# -----------------------------------------------------------------------------
num_results <- list()
for (cov in numerical_covariates) {
  vals <- suppressWarnings(as.numeric(meta_raw[[cov]]))
  if (all(is.na(vals))) {
    log_msg("  %-22s SKIPPED - no numeric values present", cov)
    next
  }
  by_study <- meta_raw %>%
    mutate(.val = vals) %>%
    filter(!is.na(.val)) %>%
    group_by(Study_ID) %>%
    summarise(n = n(), median = median(.val), q1 = quantile(.val, 0.25),
              q3 = quantile(.val, 0.75), .groups = "drop")
  out_path <- file.path(R_OUT_DIR, "confounding", sprintf("distribution_%s.tsv", cov))
  fwrite(by_study, out_path, sep = "\t")

  kw <- tryCatch(kruskal.test(vals ~ meta_raw$Study_ID), error = function(e) NULL)
  num_results[[cov]] <- list(covariate = cov, p_value = if (!is.null(kw)) kw$p.value else NA)
  log_msg("  %-22s Kruskal-Wallis p=%s across studies",
          cov, if (!is.null(kw)) sprintf("%.4f", kw$p.value) else "NA")
}

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
summary_lines <- c(
  "=== Batch confounding check summary ===",
  "",
  "Categorical covariates (chi-squared test of independence vs Study_ID):",
  sapply(cat_results, function(r) sprintf("  %-22s p=%s", r$covariate,
         ifelse(is.na(r$p_value), "NA", sprintf("%.4f", r$p_value)))),
  "",
  "Numerical covariates (Kruskal-Wallis across studies):",
  sapply(num_results, function(r) sprintf("  %-22s p=%s", r$covariate,
         ifelse(is.na(r$p_value), "NA", sprintf("%.4f", r$p_value)))),
  "",
  "A low p-value here is EXPECTED and does not by itself mean batch",
  "correction is unsafe: with many studies, almost any covariate will show",
  "some association with Study_ID by chance or by design. What matters is",
  "reading the contingency/distribution tables directly: does one or two",
  "studies dominate a covariate level, or is the association spread evenly",
  "across many studies? The former is a real confound worth discussing",
  "before batch correction; the latter is expected heterogeneity."
)
out_summary <- file.path(R_OUT_DIR, "confounding", "confounding_check_summary.txt")
writeLines(summary_lines, out_summary)
log_msg("[4/4] Summary written -> %s", out_summary)
log_msg("=== DONE ===")
