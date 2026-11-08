#!/usr/bin/env Rscript
# =============================================================================
# inspect_bactotraits.R
#
# Reports the structure of whatever BactoTraits file has been downloaded, and
# checks how many of this project's 84 genera it actually covers.
#
# RUN:  Rscript inspect_bactotraits.R /path/to/bactotraits_file
#
# -----------------------------------------------------------------------------
# WHY COVERAGE IS THE WHOLE QUESTION
# -----------------------------------------------------------------------------
# BactoTraits records MEASURED phenotypes - oxygen tolerance, Gram stain, spore
# formation, cell shape, motility - taken from BacDive, which catalogues
# CULTURED strains. Nothing is inferred from a phylogeny, which is precisely
# what makes it different from PICRUSt2 and why the failed PICRUSt2 validation
# does not rule it out.
#
# The cost of that is coverage. A genus with no cultured representative has no
# traits, and chicken caecal microbiota is heavily uncultured - the 30% of ASVs
# that could not be classified to genus in this project is a symptom of the
# same problem.
#
# If most genera match, the trait analysis is worth building. If few do, the
# analysis is not possible and that is a complete answer in itself: it can be
# reported as attempted and infeasible rather than left as an open question.
#
# This script does not assume a file layout. Trait databases vary, so it
# reports what is actually there before anything is built against it.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript inspect_bactotraits.R /path/to/file")
f <- args[1]
if (!file.exists(f)) stop("Not found: ", f)

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  ff <- grep("^--file=", a, value = TRUE)
  if (length(ff) > 0) dirname(normalizePath(sub("^--file=", "", ff))) else normalizePath(".")
}
suppressWarnings(try(source(file.path(.locate(), "00_r_config.R")), silent = TRUE))

cat("=== FILE ===\n")
cat("  ", f, "\n")
cat("   size:", round(file.info(f)$size/1024/1024, 1), "MB\n\n")

ext <- tolower(tools::file_ext(f))
d <- NULL
if (ext %in% c("xlsx","xls")) {
  sh <- readxl::excel_sheets(f)
  cat("=== SHEETS ===\n"); for (s in sh) cat("  ", s, "\n"); cat("\n")
  d <- as.data.frame(readxl::read_excel(f, sheet = 1), stringsAsFactors = FALSE)
} else {
  d <- as.data.frame(data.table::fread(f), stringsAsFactors = FALSE)
}

cat("=== DIMENSIONS ===\n")
cat("  ", nrow(d), "rows x", ncol(d), "columns\n\n")

cat("=== COLUMNS ===\n")
for (i in seq_along(names(d))) {
  v <- d[[i]]
  nu <- length(unique(v[!is.na(v)]))
  ex <- paste(head(unique(v[!is.na(v)]), 3), collapse = " | ")
  cat(sprintf("  %2d  %-34s %6d unique   e.g. %s\n", i, names(d)[i], nu, substr(ex, 1, 44)))
}

cat("\n=== LIKELY TAXONOMY COLUMNS ===\n")
taxcols <- grep("genus|taxon|species|strain|name|phylum|family|order|class",
                names(d), ignore.case = TRUE, value = TRUE)
if (length(taxcols)) for (tc in taxcols) {
  cat(sprintf("  %-30s  %d unique values\n", tc, length(unique(d[[tc]]))))
  cat("      examples:", paste(head(unique(d[[tc]][!is.na(d[[tc]])]), 5), collapse=", "), "\n")
} else cat("  none obviously named - check the column list above\n")

# -----------------------------------------------------------------------------
# coverage against this project's genera
# -----------------------------------------------------------------------------
clr <- "/scratch/user/$USER/pdi/r_analysis/genus/genus_clr_batch_adjusted.tsv"
if (file.exists(clr)) {
  g <- data.table::fread(clr)[[1]]
  gclean <- gsub("^g__", "", g)
  cat("\n=== COVERAGE AGAINST THIS PROJECT'S", length(g), "GENERA ===\n")
  best <- NULL; bestn <- -1
  for (tc in taxcols) {
    vals <- unique(trimws(as.character(d[[tc]])))
    # genus is usually the first word of a strain name
    firstword <- unique(sub("^([A-Za-z]+).*$", "\\1", vals))
    n <- sum(gclean %in% firstword)
    cat(sprintf("  matching on %-28s %3d of %d genera (%.0f%%)\n",
                tc, n, length(gclean), 100*n/length(gclean)))
    if (n > bestn) { bestn <- n; best <- tc }
  }
  cat("\n  best column:", best, "-", bestn, "of", length(gclean), "genera matched\n\n")
  if (bestn >= 0.5*length(gclean)) {
    cat("  READING: coverage is adequate. A trait analysis is worth building.\n")
  } else if (bestn >= 0.25*length(gclean)) {
    cat("  READING: partial coverage. An analysis is possible but rests on a\n")
    cat("  minority of the community, and the unmatched genera will not be a\n")
    cat("  random subset - they will be the poorly characterised ones.\n")
  } else {
    cat("  READING: coverage is too low. The trait approach is not feasible for\n")
    cat("  this corpus, which is itself a reportable result.\n")
  }
  if (bestn > 0 && bestn < length(gclean)) {
    vals <- unique(trimws(as.character(d[[best]])))
    fw <- unique(sub("^([A-Za-z]+).*$", "\\1", vals))
    miss <- setdiff(gclean, fw)
    cat("\n  unmatched genera (", length(miss), "):\n", sep="")
    cat("   ", paste(head(miss, 40), collapse=", "), "\n")
  }
} else {
  cat("\n  (genus table not found - run this on the cluster to get coverage)\n")
}
