#!/usr/bin/env Rscript
# =============================================================================
# 01_genus_agglomeration.R
#
# Collapses the ASV-level feature table to genus level using the SILVA
# taxonomy assignments, and applies a prevalence filter.
#
# HOW TO RUN: from a terminal (Rscript 01_genus_agglomeration.R) or inside
# RStudio (open the file, click Source). No path editing required; see
# 00_r_config.R for the expected repository layout.
#
# WHY GENUS LEVEL: ASVs are the correct unit for exact-sequence comparison
# within the upstream QIIME 2 pipeline, but the index itself is built at
# genus level. Two reasons: species-level SILVA assignments are documented
# as unreliable and are not used here; and genus-level agglomeration is more
# robust to the amplicon-region heterogeneity in a multi-study corpus (e.g.
# V3-V4, V4, V3-V5, full-length reads), since ASV-level differences between
# regions collapse once summed within a shared genus.
#
# PREVALENCE FILTER: a genus must be present (non-zero) in at least the
# specified fraction of samples to be retained. This removes genera observed
# in only a handful of samples, which cannot support a stable coefficient in
# a downstream regression regardless of their total abundance.
# =============================================================================

# ---- Locate this script's folder and load shared configuration -------------
# Duplicated (rather than imported) from 00_r_config.R, because that file has
# not been sourced yet at this point - this snippet's only job is to find it.
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

PREVALENCE_MIN <- 0.10   # genus must appear in >=10% of samples

# -----------------------------------------------------------------------------
# 16S REGION RESTRICTION - APPLIED HERE, NOT DOWNSTREAM
#
# The corpus spans four amplicon regions: 44 studies V3-V4 (2,784 runs), 4 V4,
# 1 V3-V5 (Lucke_2018), 1 PacBio full-length (Tian_2026) - 253 runs, 8.3% of
# the total, outside V3-V4. Primer choice determines which genera amplify at
# all, so pooling regions mixes a detection artefact into a biological signal.
# MMUPHin cannot remove it (it models a systematic shift between batches, not
# per-taxon coverage), and closed-reference clustering cannot either (it
# harmonises NAMES; it cannot recover reads never generated).
#
# WHY THIS FILTER SITS AT AGGLOMERATION RATHER THAN AT THE MODELLING STAGE
#
# Merging and taxonomy are unaffected by region: ASVs are exact sequences, so
# a V4 ASV and a V3-V4 ASV are different features and never collide. The
# mixing begins HERE, at genus collapse, where a V4 Lactobacillus ASV and a
# V3-V4 Lactobacillus ASV land in the same row. Everything downstream of that
# then inherits the mixture:
#
#   - the prevalence filter below would be computed across all 3,037 runs, so
#     a genus amplifying well under V4 but poorly under V3-V4 could clear 10%
#     on the strength of 253 foreign runs and then be near-zero in the actual
#     analysis set
#   - cmultRepl zero-replacement in 03 operates on the whole matrix
#   - MMUPHin in 03 fits its correction including those studies
#
# Filtering at the regression stage instead would remove those runs from the
# model while leaving them in the feature space, the prevalence filter, the
# zero-replacement and the batch model - a partial restriction presented as a
# clean one. Applying it here makes the restriction actually true.
#
# All 50 studies remain in the systematic review, the risk-of-bias table and
# the narrative synthesis, and were processed identically through denoising,
# merging and taxonomy. This affects index derivation only.
#
# Set REGION_FILTER to NA for the pooled-region sensitivity analysis.
# -----------------------------------------------------------------------------
REGION_FILTER <- "V3-V4"

# -----------------------------------------------------------------------------
# PRODUCTION TYPE
#
# The registered population is domestic chickens of any production type, so
# broilers, layers and breeders are all eligible and all are included by
# default. Of the 50 studies, 11 are layers and 1 is a breeder flock.
#
# They are not interchangeable birds. A layer sampled at 500 days and a broiler
# sampled at 42 differ in growth rate, gut morphology and microbiome maturity,
# and that difference is a plausible contributor to the high between-study
# heterogeneity observed in the GHI validation.
#
# Set BIRD_FILTER to "Broilers" to repeat the derivation on broiler studies
# only, as a sensitivity analysis. NA keeps every production type, which is
# the registered and primary analysis.
# -----------------------------------------------------------------------------
BIRD_FILTER <- { e <- Sys.getenv("BIRD_FILTER", ""); if (nzchar(e)) e else NA_character_ }

# -----------------------------------------------------------------------------
# STEP 1 - LOAD THE FEATURE TABLE
# The exported TSV has a comment line (starting with '#') before the header,
# which QIIME 2's biom-to-tsv export always includes; it is skipped here.
# -----------------------------------------------------------------------------
log_msg("[1/5] Loading feature table from %s", FEATURE_TABLE_TSV)
ft_raw <- fread(FEATURE_TABLE_TSV, skip = 1, header = TRUE, sep = "\t",
                check.names = FALSE, na.strings = c("", "NA"))

# First column is the feature (ASV) ID; taxonomy was attached to the last
# columns during export (upstream pipeline stage 04). Confidence, if
# present, belongs to the taxonomy assignment rather than the abundance
# matrix and is excluded from the sample columns.
feature_col <- names(ft_raw)[1]
tax_col     <- if ("taxonomy" %in% names(ft_raw)) "taxonomy" else names(ft_raw)[ncol(ft_raw)]
sample_cols <- setdiff(names(ft_raw), c(feature_col, tax_col, "confidence"))

log_msg("  %d features x %d samples", nrow(ft_raw), length(sample_cols))

# -----------------------------------------------------------------------------
# STEP 2 - PARSE TAXONOMY DOWN TO GENUS
# SILVA taxonomy strings follow: d__;p__;c__;o__;f__;g__;s__
# An ASV unclassified at genus level is labelled explicitly rather than
# silently dropped or merged into an unrelated genus - "g__unclassified" is
# a real, trackable category in the resulting table, not a loss of data.
# -----------------------------------------------------------------------------
log_msg("[2/5] Parsing taxonomy to genus level")

parse_genus <- function(tax_string) {
  parts <- str_split(tax_string, ";")[[1]] %>% str_trim()
  genus_part <- parts[str_starts(parts, "g__")]
  if (length(genus_part) == 0 || genus_part == "g__" || is.na(genus_part)) {
    return("g__unclassified")
  }
  genus_part
}

ft_raw[, genus := sapply(get(tax_col), parse_genus)]

n_unclassified <- sum(ft_raw$genus == "g__unclassified")
log_msg("  %d of %d ASVs (%.1f%%) unclassified at genus level",
        n_unclassified, nrow(ft_raw), 100 * n_unclassified / nrow(ft_raw))

# -----------------------------------------------------------------------------
# STEP 3 - AGGLOMERATE
# Sums ASV counts within each genus, per sample. This is a straightforward
# sum because feature counts are additive - unlike relative abundances,
# which would need renormalisation after collapsing.
# -----------------------------------------------------------------------------
log_msg("[3/5] Agglomerating to genus level")

# ---- region restriction, before agglomeration and before prevalence ---------
if (!is.na(REGION_FILTER)) {
  md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  md <- md[!startsWith(as.character(md[[1]]), "#")]      # drop the q2:types row
  if (!"Region_16S" %in% names(md))
    stop("Region_16S column absent from ", METADATA_TSV)

  run2reg <- setNames(trimws(as.character(md$Region_16S)),
                      trimws(as.character(md[[1]])))
  run2std <- setNames(trimws(as.character(md$Study_ID)),
                      trimws(as.character(md[[1]])))

  col_reg <- run2reg[sample_cols]
  unknown <- sample_cols[is.na(col_reg)]
  if (length(unknown))
    stop(length(unknown), " sample column(s) have no Region_16S in ",
         basename(METADATA_TSV), ", e.g. ",
         paste(utils::head(unknown, 5), collapse = ", "),
         "\n  Fix the metadata rather than letting unlabelled runs through.")

  keep_cols <- sample_cols[col_reg == REGION_FILTER]
  drop_cols <- setdiff(sample_cols, keep_cols)
  dropped_studies <- sort(unique(run2std[drop_cols]))

  log_msg("  region restriction: %s only", REGION_FILTER)
  log_msg("    %d of %d runs retained (%.1f%%)", length(keep_cols),
          length(sample_cols), 100 * length(keep_cols) / length(sample_cols))
  log_msg("    %d stud%s excluded from index derivation:", length(dropped_studies),
          if (length(dropped_studies) == 1) "y" else "ies")
  for (s in dropped_studies)
    log_msg("      %-20s %-28s %d runs", s,
            unique(run2reg[names(run2std)[run2std == s]])[1],
            sum(run2std[drop_cols] == s))

  ft_raw <- ft_raw[, c(setdiff(names(ft_raw), sample_cols), keep_cols), with = FALSE]
  sample_cols <- keep_cols
  log_msg("    prevalence, CLR and batch correction will now see %s only",
          REGION_FILTER)
} else {
  log_msg("  REGION_FILTER is NA - all regions pooled (sensitivity analysis)")
}

# ---- production type restriction (see BIRD_FILTER above) --------------------
if (!is.na(BIRD_FILTER)) {
  md2 <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  md2 <- md2[!startsWith(as.character(md2[[1]]), "#")]
  if (!"Bird_Type" %in% names(md2)) stop("Bird_Type column absent from ", METADATA_TSV)
  r2b <- setNames(trimws(as.character(md2$Bird_Type)), trimws(as.character(md2[[1]])))
  r2s <- setNames(trimws(as.character(md2$Study_ID)),  trimws(as.character(md2[[1]])))
  bt  <- r2b[sample_cols]
  keep2 <- !is.na(bt) & bt == BIRD_FILTER
  dropped_s <- sort(unique(r2s[sample_cols[!keep2]]))
  log_msg("  production type restriction: %s only", BIRD_FILTER)
  log_msg("    %d of %d runs retained; %d stud%s excluded",
          sum(keep2), length(keep2), length(dropped_s),
          if (length(dropped_s) == 1) "y" else "ies")
  for (s in dropped_s) log_msg("      %-24s %s", s, r2b[names(r2s)[r2s == s]][1])
  ft_raw <- ft_raw[, c(setdiff(names(ft_raw), sample_cols),
                       sample_cols[keep2]), with = FALSE]
  sample_cols <- sample_cols[keep2]
  log_msg("  %d runs across %d studies after restriction",
          length(sample_cols), length(unique(r2s[sample_cols])))
} else {
  log_msg("  BIRD_FILTER is NA - all production types included (primary analysis)")
}

genus_table <- ft_raw[, lapply(.SD, sum, na.rm = TRUE), by = genus, .SDcols = sample_cols]
log_msg("  %d genera before prevalence filtering", nrow(genus_table))

# -----------------------------------------------------------------------------
# STEP 4 - PREVALENCE FILTER
# -----------------------------------------------------------------------------
log_msg("[4/5] Applying prevalence filter (>= %.0f%% of samples)", 100 * PREVALENCE_MIN)

n_samples <- length(sample_cols)
prevalence <- rowSums(genus_table[, ..sample_cols] > 0) / n_samples
genus_table_filt <- genus_table[prevalence >= PREVALENCE_MIN]

log_msg("  %d genera retained (%.1f%% of pre-filter count)",
        nrow(genus_table_filt), 100 * nrow(genus_table_filt) / nrow(genus_table))

# -----------------------------------------------------------------------------
# STEP 5 - EXPORT
# -----------------------------------------------------------------------------
log_msg("[5/5] Writing output")

out_counts <- file.path(R_OUT_DIR, "genus", "genus_counts.tsv")
fwrite(genus_table_filt, out_counts, sep = "\t")
log_msg("  genus count table -> %s", out_counts)

summary_stats <- list(
  n_asvs_input              = nrow(ft_raw),
  n_asvs_unclassified_genus = n_unclassified,
  n_genera_before_filter    = nrow(genus_table),
  n_genera_after_filter     = nrow(genus_table_filt),
  prevalence_threshold      = PREVALENCE_MIN,
  n_samples                 = n_samples
)
out_summary <- file.path(R_OUT_DIR, "genus", "genus_agglomeration_summary.txt")
writeLines(sprintf("%s: %s", names(summary_stats), unlist(summary_stats)), out_summary)
log_msg("  summary -> %s", out_summary)

log_msg("=== DONE ===")
