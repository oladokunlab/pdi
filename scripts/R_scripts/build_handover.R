#!/usr/bin/env Rscript
# =============================================================================
# build_handover.R
#
# Assemble a self-contained folder that lets someone else apply the GHI and
# the GDI without access to this pipeline.
#
# RUN:  Rscript build_handover.R          (needs 04b and 07 to have run)
# OUT:  <repo>/handover/
#
# The folder contains the two index definitions, the two scoring scripts, a
# README, and an example input for each index drawn from the derivation corpus.
# The examples serve to confirm the scripts execute correctly and to document
# the expected input format.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
SD <- .locate()
source(file.path(SD, "00_r_config.R"))

HO <- file.path(PDI_ROOT, "handover")
dir.create(HO, recursive = TRUE, showWarnings = FALSE)
log_msg("Building handover package in %s", HO)

# ---- definitions -------------------------------------------------------------
need <- c(ghi = file.path(R_OUT_DIR, "ghi", "ghi_index_definition.tsv"),
          gdi = file.path(R_OUT_DIR, "gdi", "gdi_index_definition.tsv"))
for (n in names(need)) {
  if (!file.exists(need[[n]]))
    stop("Missing ", need[[n]], "\n  Run ",
         if (n == "ghi") "04b_ghi_pca.R" else "07_gdi_construct.R", " first.")
  file.copy(need[[n]], HO, overwrite = TRUE)
  log_msg("  %s", basename(need[[n]]))
}

# ---- scripts -----------------------------------------------------------------
for (s in c("score_new_ghi.R", "score_new_samples.R")) {
  p <- file.path(SD, s)
  if (!file.exists(p)) stop("Missing script: ", p)
  file.copy(p, HO, overwrite = TRUE)
  log_msg("  %s", s)
}

# ---- worked example, GHI -----------------------------------------------------
# Studies are ranked by how many domains are jointly scorable in their rows,
# rather than by marker count alone. A study can hold many markers while no
# single row carries a complete domain, which produces an example where most
# rows score on one domain or none, and demonstrates nothing.
bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
defn_g <- fread(file.path(HO, "ghi_index_definition.tsv"))
mk <- intersect(unique(defn_g$marker), names(bio))

# domains scorable per row: a PCA domain needs all its markers, a mean domain
# needs at least one
n_dom <- rep(0L, nrow(bio))
for (dom in unique(defn_g$domain)) {
  sub <- defn_g[domain == dom]
  cols <- intersect(sub$marker, mk)
  if (!length(cols)) next
  M <- as.matrix(bio[, cols, drop = FALSE])
  ok <- if (sub$method[1] == "PCA" && length(cols) == nrow(sub))
          stats::complete.cases(M) else rowSums(!is.na(M)) >= 1
  n_dom <- n_dom + as.integer(ok)
}
bio$.n_dom <- n_dom

cand <- data.table(Study_ID = bio$Study_ID, n = n_dom)[
  , .(rows = sum(n >= 2), best = max(n)), by = Study_ID][order(-rows, -best)]
pick <- cand$Study_ID[seq_len(min(2, nrow(cand)))]

id_cols <- intersect(c("Study_ID", "Sample_ID", "Treatment_Group", "Segment"),
                     names(bio))
ex_ghi <- bio[bio$Study_ID %in% pick & bio$.n_dom >= 2, c(id_cols, mk), drop = FALSE]
if (nrow(ex_ghi) < 4)
  ex_ghi <- bio[bio$Study_ID %in% pick & bio$.n_dom >= 1, c(id_cols, mk), drop = FALSE]
fwrite(ex_ghi, file.path(HO, "example_ghi_input.tsv"), sep = "\t")
log_msg("  example_ghi_input.tsv  (%d rows, studies: %s, %d+ domains per row)",
        nrow(ex_ghi), paste(pick, collapse = ", "),
        min(bio$.n_dom[bio$Study_ID %in% pick & bio$.n_dom >= 1]))

# ---- worked example, GDI -----------------------------------------------------
# The eight panel genera for a small number of runs, in CLR space.
clr_p <- file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")
if (file.exists(clr_p)) {
  clr <- fread(clr_p)
  gcol <- names(clr)[1]
  panel <- unique(fread(file.path(HO, "gdi_index_definition.tsv"))$genus)
  have <- intersect(panel, clr[[gcol]])
  if (length(have) == length(panel)) {
    cols <- names(clr)[2:min(7, ncol(clr))]
    ex_gdi <- clr[match(panel, clr[[gcol]]), c(gcol, cols), with = FALSE]
    fwrite(ex_gdi, file.path(HO, "example_gdi_input.tsv"), sep = "\t")
    log_msg("  example_gdi_input.tsv  (%d genera x %d samples)",
            nrow(ex_gdi), length(cols))
  } else {
    log_msg("  example_gdi_input.tsv SKIPPED - %d panel genera absent from the CLR table",
            length(setdiff(panel, have)))
  }
} else {
  log_msg("  example_gdi_input.tsv SKIPPED - %s not found", basename(clr_p))
}

# ---- README ------------------------------------------------------------------
rd <- file.path(SD, "..", "handover_README.md")
if (file.exists(rd)) {
  file.copy(rd, file.path(HO, "README.md"), overwrite = TRUE)
  log_msg("  README.md")
} else {
  log_msg("  NOTE: handover_README.md not found beside the scripts folder;")
  log_msg("        copy it into %s manually.", HO)
}

log_msg("")
log_msg("Verify with:")
log_msg("  cd %s", HO)
log_msg("  Rscript score_new_ghi.R example_ghi_input.tsv")
log_msg("  Rscript score_new_samples.R example_gdi_input.tsv")
log_msg("=== DONE ===")
