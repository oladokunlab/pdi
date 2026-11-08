#!/bin/bash
#SBATCH --job-name=pdi_merge_tax
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --output=/scratch/user/%u/pdi/logs/merge_tax_%j.log
# =============================================================================
# 03_merge_and_taxonomy.sh - merge all studies, build tree, assign taxonomy
#
# HOW TO RUN: sbatch, only after stage 02 has finished for every study.
#   sbatch 03_merge_and_taxonomy.sh
#
# MEMORY WARNING: classify-sklearn against the full-length SILVA classifier is
# genuinely memory-hungry at corpus scale. 180G is allocated deliberately;
# if this still runs out of memory, reduce --p-reads-per-batch rather than
# increasing memory further.
# =============================================================================

# NOTE: under sbatch, SLURM copies this script to /var/spool/slurmd/jobNNNN/
# and executes it there, so $(dirname "$0") points at the spool directory, NOT
# at the scripts directory, and sourcing 00_config.sh would fail.
# SLURM_SUBMIT_DIR is the directory sbatch was invoked from, which is where
# the repository actually lives.
# Falls back to $0's directory when run interactively on the login node.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
if [[ ! -f "${SCRIPT_DIR}/00_config.sh" ]]; then
  echo "FATAL: 00_config.sh not found in ${SCRIPT_DIR}" >&2
  echo "Submit with:  cd /scratch/user/$USER/pdi/scripts && sbatch $(basename "$0")" >&2
  exit 1
fi
source "${SCRIPT_DIR}/00_config.sh"

# hard stop if config did not actually load - otherwise every downstream
# variable is empty and the job burns SUs doing nothing
if [[ -z "${PDI_ROOT:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh loaded but PDI_ROOT/log() are unset." >&2
  exit 1
fi
activate_qiime

cd "$MERGED_DIR"

# =============================================================================
# STEP 1 - MERGE FEATURE TABLES AND REP-SEQS
# This is where the separately-denoised studies become one dataset. Valid
# because ASVs are exact sequences - identical sequence = identical feature ID
# (an md5 hash of the sequence), so the merge lines up correctly across studies.
# Samples present in only one study simply get zeros everywhere else.
# =============================================================================
if ! already_done "${MERGED_DIR}/merged-table.qza"; then
  log "[1/4] Merging per-study feature tables"

  TABLES=(); SEQS=()
  while IFS=$'\t' read -r study n region platform layout pf pr status; do
    [[ "$status" == HOLD_* ]] && continue
    t="${WORK_DIR}/${study}/table.qza"
    s="${WORK_DIR}/${study}/rep-seqs.qza"
    if [[ -s "$t" && -s "$s" ]]; then
      TABLES+=(--i-tables "$t"); SEQS+=(--i-data "$s")
    else
      log "WARN: ${study} has no denoised output - EXCLUDED from merge."
    fi
  done < <(tail -n +2 "$STUDY_TABLE")

  log "Merging $(( ${#TABLES[@]} / 2 )) study tables"

  $QIIME feature-table merge "${TABLES[@]}" \
    --p-overlap-method error_on_overlapping_sample \
    --o-merged-table "${MERGED_DIR}/merged-table.qza"
  # error_on_overlapping_sample is intentional: if two studies somehow contain
  # the same sample ID, this produces a hard failure rather than a silent
  # overwrite that would corrupt the merged table.

  $QIIME feature-table merge-seqs "${SEQS[@]}" \
    --o-merged-data "${MERGED_DIR}/merged-rep-seqs.qza"
fi

# =============================================================================
# STEP 2 - SUMMARISE THE MERGED TABLE
# This .qzv is where stage 04 reads the MEAN SAMPLE DEPTH from, for the
# rare-ASV threshold in stage 04. It also provides an initial view of the
# corpus as a whole:
# how many samples, how many ASVs, how uneven the depth is.
# =============================================================================
if ! already_done "${MERGED_DIR}/merged-table-summary.qzv"; then
  log "[2/4] Summarising merged table"
  $QIIME feature-table summarize \
    --i-table "${MERGED_DIR}/merged-table.qza" \
    --m-sample-metadata-file "$METADATA" \
    --o-visualization "${MERGED_DIR}/merged-table-summary.qzv"
fi

# =============================================================================
# STEP 3 - TAXONOMY (SILVA 138, 99%, FULL-LENGTH)
# A naive Bayes classifier trained on SILVA. It takes each ASV sequence and
# returns a taxonomic string (kingdom -> genus, sometimes species) plus a
# confidence value.
#
# FULL-LENGTH classifier, chosen deliberately: region-specific classifiers are
# trained on reads extracted between one primer pair. Mine is a mixed-region
# corpus, so no single region-trimmed classifier fits. The full-length model
# classifies all of them consistently.
#
# --p-read-orientation auto is used instead of 'same', which is appropriate for a single
# in-house run where you control library prep. Mine is 54 labs' worth of SRA
# submissions with inconsistent read orientation, so 'auto' lets QIIME2 detect
# orientation per batch. Flagged in README as a documented deviation.
# =============================================================================
if ! already_done "${MERGED_DIR}/taxonomy.qza"; then
  log "[3/4] Classifying against SILVA 138 full-length"
  [[ -s "$CLASSIFIER" ]] || { log "ERROR: classifier not at $CLASSIFIER"; exit 1; }

  $QIIME feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "${MERGED_DIR}/merged-rep-seqs.qza" \
    --p-read-orientation auto \
    --p-reads-per-batch 10000 \
    --p-n-jobs "$NCORES" \
    --o-classification "${MERGED_DIR}/taxonomy.qza"
fi

# =============================================================================
# STEP 5 - VIEWABLE TAXONOMY + EXPORT
# =============================================================================
if ! already_done "${MERGED_DIR}/taxonomy.qzv"; then
  log "[4/4] Tabulating and exporting taxonomy"
  $QIIME metadata tabulate \
    --m-input-file "${MERGED_DIR}/taxonomy.qza" \
    --o-visualization "${MERGED_DIR}/taxonomy.qzv"

  $QIIME tools export \
    --input-path "${MERGED_DIR}/taxonomy.qza" \
    --output-path "${MERGED_DIR}/taxa"
fi

log "=== STAGE 03 COMPLETE ==="
log "Next: open merged-table-summary.qzv, then run 04_filter_and_export.sh"
