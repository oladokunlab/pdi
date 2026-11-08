#!/usr/bin/env Rscript
# =============================================================================
# 00_r_config.R
#
# Shared paths and parameters for the downstream analysis scripts (genus
# agglomeration, batch confounding check, GHI/GDI construction). Sourced by
# every other script in this folder. Edit paths and parameters here only.
#
# EXPECTED REPOSITORY LAYOUT (paths below are computed relative to this
# file's own location, so the repository can be cloned or moved anywhere):
#
#   <repo root>/
#   |-- merged/
#   |   `-- feature-table_w_tax.tsv     genus-annotated feature table
#   |-- scripts/
#   |   |-- GDI_metadata_qiime2.tsv     QIIME 2 formatted sample metadata
#   |   |-- GDI_Metadata.xlsx           biomarker workbook (GHI construction)
#   |   `-- R_scripts/                  <- this script's folder
#   |       |-- 00_r_config.R
#   |       |-- 01_genus_agglomeration.R
#   |       `-- 02b_batch_confounding_check.R
#   `-- r_analysis/                     created automatically; all output
#                                        from these scripts is written here
#
# HOW TO RUN: works identically from the terminal (Rscript 01_....R) or from
# inside RStudio (open the file and click Source), on Windows, macOS, or
# Linux, without editing any path by hand.
# =============================================================================

# ---- Locate this script's own folder, regardless of how it was launched ----
# Three execution contexts are handled, in order of likelihood:
#   1. Rscript from a terminal      -> the file path is in commandArgs()
#   2. RStudio "Source" button      -> rstudioapi reports the open file's path
#   3. Anything else (e.g. the script was pasted into a live R console)
#      -> falls back to the current working directory, with a warning, since
#         no method exists in base R alone to determine the file's location.
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path)))
    }
  }
  warning(
    "Could not auto-detect this script's folder (not run via Rscript or ",
    "RStudio's Source button). Falling back to the current working ",
    "directory. If path errors follow, setwd() to the R_scripts folder ",
    "before running, or open this file in RStudio and use Source."
  )
  normalizePath(".")
}

SCRIPT_DIR <- get_script_dir()

# Repository root is two levels above R_scripts/ (R_scripts -> scripts -> root).
PDI_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

# ---- Derived paths -----------------------------------------------------------
MERGED_DIR <- file.path(PDI_ROOT, "merged")
R_OUT_DIR  <- file.path(PDI_ROOT, "r_analysis")

dir.create(R_OUT_DIR,                               showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(R_OUT_DIR, "genus"),           showWarnings = FALSE)
dir.create(file.path(R_OUT_DIR, "confounding"),     showWarnings = FALSE)
dir.create(file.path(R_OUT_DIR, "ghi"),             showWarnings = FALSE)

# ---- Input files ---------------------------------------------------------------
# The QIIME 2 pipeline writes the annotated feature table into
# merged/final_exported/ (stage 04's export step). Some setups place it
# directly in merged/. Both are checked here, in that order, so the script
# works without hand-editing this path on a new machine - and fails with a
# clear message naming both locations if neither exists, rather than a bare
# "file not found" from deep inside fread().
.candidate_ft <- c(
  file.path(MERGED_DIR, "final_exported", "feature-table_w_tax.tsv"),
  file.path(MERGED_DIR, "feature-table_w_tax.tsv")
)
FEATURE_TABLE_TSV <- .candidate_ft[file.exists(.candidate_ft)][1]
if (is.na(FEATURE_TABLE_TSV)) {
  stop(sprintf(
    "Annotated feature table not found. Looked in:\n  %s\nRun stage 04 (04_filter_and_export.sh) first.",
    paste(.candidate_ft, collapse = "\n  ")))
}
# The metadata lives alongside the scripts. SCRIPT_DIR is R_scripts/, so its
# parent is the scripts folder - derived rather than hardcoded as "scripts",
# so that directory can be renamed (e.g. GDI_pipeline) without breaking this.
SCRIPTS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
METADATA_TSV      <- file.path(SCRIPTS_DIR, "GDI_metadata_qiime2.tsv")
METADATA_XLSX     <- file.path(SCRIPTS_DIR, "GDI_Metadata.xlsx")
BIOMARKER_SHEET    <- "Sample_with_Zscores"

# ---- 16S region, from the same studies.tsv that drove denoising ---------------
# Region is a study-level property that until now lived only on the bash side,
# which is why the V3-V4 restriction could not be applied in R. Reading the
# same file the SLURM array used guarantees the two sides cannot drift apart.
STUDIES_TSV <- file.path(SCRIPTS_DIR, "studies.tsv")

study_regions <- function() {
  if (!file.exists(STUDIES_TSV)) stop("Not found: ", STUDIES_TSV)
  s <- data.table::fread(STUDIES_TSV, sep = "\t", header = TRUE)
  setNames(trimws(as.character(s$region)), trimws(as.character(s$study_id)))
}

# ---- Logging -------------------------------------------------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
}

# ---- Package check ---------------------------------------------------------------
# Lists only the packages the current scripts actually call, and checks each
# by name individually rather than requiring the full tidyverse meta-package,
# since a fresh R installation typically has none of these pre-installed.
# Fails loudly with an exact install command, rather than failing later
# inside a specific analysis step with a less informative error.
required_packages <- c("data.table", "readxl", "stringr", "dplyr", "tidyr")
missing_pkgs <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(sprintf(
    "Missing required R packages: %s\nInstall with: install.packages(c(%s))",
    paste(missing_pkgs, collapse = ", "),
    paste(sprintf('"%s"', missing_pkgs), collapse = ", ")
  ))
}

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

log_msg("R config loaded. Repository root detected at: %s", PDI_ROOT)
