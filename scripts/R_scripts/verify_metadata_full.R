#!/usr/bin/env Rscript
# =============================================================================
# verify_metadata_full.R
#
# Systematic pre-run check of BOTH metadata files. Read-only, writes nothing.
# Checks the things that actually break the pipeline downstream, not just
# whether a column looks tidy.
# =============================================================================
.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
PASS <- 0; FAIL <- 0
chk <- function(cond, m, ...) {
  msg <- sprintf(m, ...)
  if (cond) { log_msg("  PASS  %s", msg); PASS <<- PASS + 1 }
  else      { log_msg("  FAIL  %s", msg); FAIL <<- FAIL + 1 }
}

log_msg("=== BIOMARKER WORKBOOK ===")
sheets <- readxl::excel_sheets(METADATA_XLSX)
bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET), stringsAsFactors = FALSE)
bio$Age_d <- suppressWarnings(as.numeric(bio$Age_d))
log_msg("  %d sheets, biomarker sheet: %d rows x %d cols", length(sheets), nrow(bio), ncol(bio))

chk(nrow(bio) == 566, "row count is 566 (expected)")
chk(sum(blank(bio$Age_d)) == 0, "no blank Age_d (%d blank)", sum(blank(bio$Age_d)))
chk(sum(blank(bio$Study_ID)) == 0, "no blank Study_ID (%d blank)", sum(blank(bio$Study_ID)))
chk(sum(blank(bio$Treatment_Group)) == 0, "no blank Treatment_Group (%d blank)", sum(blank(bio$Treatment_Group)))
chk(sum(blank(bio$Segment)) == 0, "no blank Segment (%d blank)", sum(blank(bio$Segment)))
chk(all(bio$Age_d >= 0 & bio$Age_d <= 600, na.rm=TRUE), "all Age_d within 0-600 days")

# duplicate keys: same study+arm+segment appearing twice is what forces the
# "averaging duplicate keys" warning downstream
key <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), trimws(bio$Segment), sep="||")
ndup <- sum(duplicated(key))
log_msg("  NOTE  %d duplicate Study+Arm+Segment key(s) - these get averaged at join time", ndup)

log_msg("")
log_msg("=== QIIME2 METADATA ===")
raw <- readLines(METADATA_TSV, warn = FALSE)
hdr <- strsplit(raw[1], "\t")[[1]]
chk(hdr[1] == "sample-id", "line 1 is the real header (starts with 'sample-id')")
chk(startsWith(raw[2], "#q2:types"), "line 2 is the #q2:types declaration")
chk(!any(startsWith(raw[-2], "#q2:types")), "no duplicate #q2:types line elsewhere")
dat <- raw[-(1:2)]
lens <- lengths(strsplit(dat, "\t", fixed=TRUE))
chk(all(lens == length(hdr)), "all %d data rows have %d fields", length(dat), length(hdr))
log_msg("  %d data rows, %d columns", length(dat), length(hdr))
chk(!any(grepl("^Age", hdr)), "no Age column here (age lives only in the workbook - correct)")

sid <- which(hdr == "Study_ID")
seq_studies <- unique(trimws(vapply(strsplit(dat,"\t",fixed=TRUE), function(r) r[sid], character(1))))

log_msg("")
log_msg("=== CROSS-FILE AGREEMENT ===")
bio_studies <- unique(trimws(bio$Study_ID))
only_bio <- setdiff(bio_studies, seq_studies)
only_seq <- setdiff(seq_studies, bio_studies)
log_msg("  biomarker studies: %d, sequencing studies: %d", length(bio_studies), length(seq_studies))
chk(length(only_seq) == 0, "every sequencing study has biomarker data (%d missing)", length(only_seq))
if (length(only_seq)) for (s in only_seq) log_msg("        sequencing-only: %s", s)
log_msg("  NOTE  biomarker-only (expected: Ma_2024_Layers): %s", paste(only_bio, collapse=", "))

# do arm names agree for studies present in both?
mismatch_arms <- character(0)
for (s in intersect(bio_studies, seq_studies)) {
  b <- unique(trimws(bio$Treatment_Group[trimws(bio$Study_ID)==s]))
  q <- unique(trimws(vapply(strsplit(dat,"\t",fixed=TRUE), function(r) r[sid], character(1))))
  qi <- which(trimws(vapply(strsplit(dat,"\t",fixed=TRUE), function(r) r[sid], character(1))) == s)
  tg <- which(hdr == "Treatment_Group")
  qarms <- unique(trimws(vapply(strsplit(dat[qi],"\t",fixed=TRUE), function(r) r[tg], character(1))))
  orphan <- setdiff(qarms, b)
  if (length(orphan)) mismatch_arms <- c(mismatch_arms, sprintf("%s: %s", s, paste(orphan, collapse=", ")))
}
chk(length(mismatch_arms) == 0, "all sequencing arm names match a biomarker arm (%d studies with orphans)", length(mismatch_arms))
if (length(mismatch_arms)) for (m in head(mismatch_arms, 15)) log_msg("        %s", m)

log_msg("")
log_msg("=========================================")
log_msg("PASS: %d    FAIL: %d", PASS, FAIL)
log_msg("=========================================")
if (FAIL == 0) log_msg("Metadata is consistent. Safe to run the pipeline.")
