#!/usr/bin/env Rscript
# =============================================================================
# v2_01_genus_agglomeration.R
#
# Genus agglomeration for the v2 (QIIME2 2026.1, project-trained V3-V4
# classifier) corpus. Mirrors 01_genus_agglomeration.R exactly, pointed at
# the v2 output tree instead.
#
# RUN (as part of a batch job, not interactively - the input file is 1.6GB):
#   Rscript v2_01_genus_agglomeration.R
#
# -----------------------------------------------------------------------------
# WHAT DIFFERS FROM THE ORIGINAL SCRIPT, AND WHY
# -----------------------------------------------------------------------------
# The original script applies a region filter (REGION_FILTER = "V3-V4")
# because its input corpus mixes four amplicon regions. The v2 corpus never
# had that mixture: v2_01_denoise.sh only denoised studies whose Region_16S
# was V3-V4 in the first place, so every study reaching this file already
# qualifies. No region filtering is applied here - not because the step was
# skipped, but because the corpus was already restricted three steps
# upstream, at the point where restricting it actually matters (before
# thresholds are computed from the data).
#
# Everything else - prevalence filtering, genus-level summing from the taxon
# string - is identical to the original.
#
# The metadata file is the SAME file the original pipeline uses
# (GDI_metadata_qiime2.tsv). Study_ID, Treatment_Group and Challenge_Type
# describe the biological sample and experimental design, not which
# bioinformatics pipeline processed its reads, so no new metadata was built
# for this run.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

V2_ROOT <- "/scratch/user/$USER/pdi/v2_2026.1"
V2_FEATURE_TABLE <- file.path(V2_ROOT, "merged/final_exported/feature-table_w_tax.tsv")
V2_OUT_DIR <- file.path(V2_ROOT, "r_analysis/genus")
dir.create(V2_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREVALENCE_MIN <- 0.10   # identical to the original script

log_msg("=== v2 genus agglomeration ===")
log_msg("Input : %s", V2_FEATURE_TABLE)
log_msg("Output: %s", V2_OUT_DIR)

if (!file.exists(V2_FEATURE_TABLE)) {
  log_msg("FATAL: input file not found at %s", V2_FEATURE_TABLE)
  quit(status = 1)
}

log_msg("[1/5] Reading feature table (1.6GB - this is the slow step)")
# Line 1 is "# Constructed from biom file", line 2 is the real header.
raw <- fread(V2_FEATURE_TABLE, skip = 1, header = TRUE, sep = "\t",
             check.names = FALSE, showProgress = TRUE)
log_msg("  loaded: %d rows x %d columns", nrow(raw), ncol(raw))

id_col <- names(raw)[1]
tax_col <- names(raw)[ncol(raw)]
log_msg("  ID column: %s | taxonomy column: %s", id_col, tax_col)
if (!grepl("tax", tax_col, ignore.case = TRUE)) {
  log_msg("WARNING: last column '%s' does not look like a taxonomy column by name.", tax_col)
  log_msg("  First few values: %s", paste(head(raw[[tax_col]], 3), collapse = " | "))
  log_msg("  Proceeding on the assumption it is taxonomy, since it is the documented")
  log_msg("  biom-export layout, but this is worth a human eye before trusting the output.")
}

sample_cols <- setdiff(names(raw), c(id_col, tax_col))
log_msg("  %d sample columns detected", length(sample_cols))

log_msg("[2/5] Extracting genus from the taxonomy string (g__... token)")
genus <- sub(".*g__([^;]*).*", "\\1", raw[[tax_col]])
genus <- trimws(genus)

# "Incertae_Sedis" is a formal taxonomic placeholder (Latin: "of uncertain
# placement"). The classifier could not resolve these to a genus and used
# this label instead of leaving the field blank. Checked directly against
# the actual data: 42,791 of 148,535 features (28.8%) carry this label,
# making it by far the single most common value in the genus slot - larger
# than the next most frequent genus by more than 3x. Counting it as a
# genuine genus assignment would understate the true unclassified rate by
# roughly 29 percentage points, which is exactly what happened on the first
# run of this script and produced a 0.0% figure that did not survive a
# direct check of the underlying taxonomy strings.
UNRESOLVED_LABELS <- c("", "Incertae_Sedis", "uncultured", "unclassified", "unidentified")
genus[genus %in% UNRESOLVED_LABELS | is.na(genus)] <- "Unclassified"
genus <- ifelse(genus == "Unclassified", "g__Unclassified", paste0("g__", genus))
n_unclassified <- sum(genus == "g__Unclassified")
log_msg("  %d of %d features (%.1f%%) unclassified or unresolved at genus level",
        n_unclassified, length(genus), 100 * n_unclassified / length(genus))
log_msg("  (includes blank, Incertae_Sedis, uncultured, unclassified, unidentified -")
log_msg("  labels the classifier uses when it cannot resolve a confident genus call)")

log_msg("[3/5] Summing counts within each genus")
mat <- as.matrix(raw[, ..sample_cols])
genus_sum <- rowsum(mat, group = genus)
log_msg("  %d distinct genera before prevalence filtering", nrow(genus_sum))

log_msg("[4/5] Prevalence filter (present in >= %.0f%% of samples)", PREVALENCE_MIN * 100)
prevalence <- rowSums(genus_sum > 0) / ncol(genus_sum)
keep <- prevalence >= PREVALENCE_MIN
genus_final <- genus_sum[keep, , drop = FALSE]
log_msg("  %d -> %d genera after prevalence filtering", nrow(genus_sum), nrow(genus_final))

log_msg("[5/5] Writing output")
out <- data.frame(Genus = rownames(genus_final), genus_final, check.names = FALSE)
fwrite(out, file.path(V2_OUT_DIR, "genus_counts.tsv"), sep = "\t")

log_msg("")
log_msg("=== DONE ===")
log_msg("Genera: %d | Samples: %d | Unclassified rate: %.1f%%",
        nrow(genus_final), ncol(genus_final), 100 * n_unclassified / length(genus))
log_msg("Output -> %s", file.path(V2_OUT_DIR, "genus_counts.tsv"))
log_msg("Next: v2_03_batch_correction.R")
