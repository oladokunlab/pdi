#!/bin/bash
# =============================================================================
# 00_config.sh
#
# Shared configuration for the GDI/PDI 16S rRNA processing pipeline.
# Sourced by every other script in this repository. Edit paths and
# parameters here only; do not hardcode them elsewhere.
#
# Project : Poultry Gut Dysbiosis Index (GDI) / Poultry Dysbiosis Index (PDI)
# Registry: PROSPERO CRD420261384137
# Cluster : Grace HPRC (Texas A&M)
# =============================================================================

set -euo pipefail

# ---- Locate this config file, so nothing depends on the folder's NAME ------
# BASH_SOURCE[0] is this file even when sourced from another script. Deriving
# the scripts directory this way means the folder can be renamed (e.g. from
# "scripts" to "GDI_pipeline") without breaking any path below.
_CFG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPTS_DIR="$_CFG_DIR"

# ---- Root paths -------------------------------------------------------------
export PDI_ROOT="/scratch/user/$USER/pdi"
export RAW_DIR="${PDI_ROOT}/raw_data"          # <Study_ID>/<RUN>_1.fastq.gz, etc.
export MANIFEST_DIR="${PDI_ROOT}/manifests"
export WORK_DIR="${PDI_ROOT}/qiime2_final"     # per-study denoising output
export MERGED_DIR="${PDI_ROOT}/merged"         # corpus-level output
export LOG_DIR="${PDI_ROOT}/logs"
export QC_DIR="${PDI_ROOT}/fastqc_out"

# ---- Metadata ---------------------------------------------------------------
# QIIME 2-formatted sample metadata. Row 1 is the header, row 2 is the
# #q2:types declaration, remaining rows are one per sequencing run.
# Checked next to the scripts first (where the distributed package puts it),
# then the older metadata/ location, so either layout works.
if [[ -s "${_CFG_DIR}/GDI_metadata_qiime2.tsv" ]]; then
  export METADATA="${_CFG_DIR}/GDI_metadata_qiime2.tsv"
else
  export METADATA="${PDI_ROOT}/metadata/GDI_metadata_qiime2.tsv"
fi

# Per-study configuration table: study_id, n_runs, region, platform, layout,
# primer_f, primer_r, status, primers_present, read_orientation. Sorted
# smallest-study-first so array jobs surface failures early rather than after
# the largest study has already run.
export STUDY_TABLE="${SCRIPTS_DIR}/studies.tsv"

# ---- QIIME 2 environment ----------------------------------------------------
# Pinned to QIIME 2 2024.10.1. The SILVA classifier below was trained with
# scikit-learn 1.4.2, which matches this build. Upgrading QIIME 2 without
# retraining the classifier will cause classify-sklearn to fail on a
# scikit-learn version mismatch.
export QIIME2_ENV="/scratch/user/$USER/.conda/envs/qiime2"
export QIIME="${QIIME2_ENV}/bin/qiime"

# SILVA 138, 99% OTUs, full-length naive Bayes classifier. A full-length
# classifier is used rather than a region-trimmed one because the corpus
# spans multiple 16S regions (V3-V4, V4, V3-V5) plus PacBio full-length
# reads; a region-trimmed classifier would be invalid for the other regions.
export CLASSIFIER="${PDI_ROOT}/classifier/silva-138-99-nb-classifier.qza"

# ---- Denoising parameters (registered in PROSPERO; do not change without ---
# ---- PI sign-off) ------------------------------------------------------------
# --p-trunc-len 0: no fixed-length truncation. Read lengths vary across the
# corpus (different platforms, different regions), so a single fixed cutoff
# is not appropriate across all studies. Quality control is instead handled
# by the expected-error model (--p-max-ee) at denoising time.
export TRUNC_LEN=0
export TRIM_LEFT=0
export MAX_EE=2.0
export CHIMERA_METHOD="consensus"

# ---- Rare-ASV filtering -----------------------------------------------------
# Minimum-frequency threshold = mean per-sample sequencing depth x this
# multiplier, rounded to the nearest integer. Computed at stage 04 from the
# merged feature table.
export RARE_ASV_MULTIPLIER=0.001

# ---- Compute -----------------------------------------------------------------
export NCORES="${SLURM_CPUS_PER_TASK:-16}"

# ---- Helpers ------------------------------------------------------------------
mkdir -p "$MANIFEST_DIR" "$WORK_DIR" "$MERGED_DIR" "$LOG_DIR" "$QC_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# Every stage checks for its own expected output before doing work, making
# the pipeline resumable: an interrupted job re-processes only the units
# that have not yet completed.
already_done() {
  [[ -s "$1" ]] && { log "SKIP (exists): $1"; return 0; }
  return 1
}

# Activates the QIIME 2 environment for the current shell. Module loading is
# attempted but not required - if no matching module is available, the
# environment's bin directory is placed on PATH directly, since the qiime
# launcher script does not depend on the module being loaded.
activate_qiime() {
  module load Miniconda3/23.10.0-1 2>/dev/null \
    || module load Anaconda3 2>/dev/null \
    || true
  export PATH="${QIIME2_ENV}/bin:$PATH"
  [[ -x "$QIIME" ]] || { echo "ERROR: qiime not found at $QIIME" >&2; exit 1; }
  local v
  v=$("$QIIME" --version 2>&1 | head -1) || {
    echo "ERROR: qiime found but failed to run: $v" >&2; exit 1; }
  log "QIIME2 OK: $v"
}
