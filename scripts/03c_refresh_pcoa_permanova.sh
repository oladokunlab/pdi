#!/bin/bash
#SBATCH --job-name=pdi_pcoa_refresh
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=/scratch/user/%u/pdi/logs/pcoa_refresh_%j.log
# =============================================================================
# 03c_refresh_pcoa_permanova.sh
#
# Regenerates ONLY the metadata-dependent QIIME2 outputs after the metadata
# corrections (Ma_2024 split into two studies, arm-name standardisation,
# age corrections).
#
# HOW TO RUN:
#   cd /scratch/user/$USER/pdi/scripts && sbatch 03c_refresh_pcoa_permanova.sh
#
# -----------------------------------------------------------------------------
# WHAT IS *NOT* RECOMPUTED, AND WHY THAT IS CORRECT
# -----------------------------------------------------------------------------
# The following are derived from SEQUENCES and COUNTS, not from metadata, and
# nothing in the metadata corrections touched either. Recomputing them would
# burn hours of compute to produce byte-identical files:
#
#   merged-table.qza, merged-rep-seqs.qza   (the merge itself)
#   rarefied-table.qza                      (rarefaction of that table)
#   braycurtis_dm.qza, jaccard_dm.qza       (distances between samples)
#   braycurtis_pcoa.qza, jaccard_pcoa.qza   (the ordination coordinates)
#   aligned-rep-seqs.qza, masked-aligned-rep-seqs.qza, rooted-tree.qza
#                                           (the phylogeny - built from
#                                            sequence alignment only)
#   taxonomy.qza                            (classifier output per ASV)
#
# WHAT *IS* REGENERATED, because it reads the metadata file:
#
#   *_emperor.qzv          the PCoA plots - points are the same, but the
#                          colour-by fields (Study_ID, Intervention_Class,
#                          PCoA_Group) come from metadata
#   *_permanova.qzv        the statistical tests - these GROUP samples by a
#                          metadata column, so the Study_ID test in particular
#                          changes now that Ma_2024 is two separate studies
#   merged-table-summary.qzv  sample metadata shown alongside the summary
#
# The existing versions of these are moved aside rather than deleted, so the
# pre-correction versions remain available for comparison.
# =============================================================================

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
if [[ ! -f "${SCRIPT_DIR}/00_config.sh" ]]; then
  echo "FATAL: 00_config.sh not found in ${SCRIPT_DIR}" >&2
  echo "Submit with:  cd /scratch/user/$USER/pdi/scripts && sbatch $(basename "$0")" >&2
  exit 1
fi
source "${SCRIPT_DIR}/00_config.sh"
if [[ -z "${PDI_ROOT:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: config loaded but PDI_ROOT/log() unset." >&2; exit 1
fi
activate_qiime

cd "$MERGED_DIR"

STAMP=$(date +%Y%m%d_%H%M)
OLD="${MERGED_DIR}/pre_metadata_correction_${STAMP}"
mkdir -p "$OLD"

log "=== REFRESHING METADATA-DEPENDENT OUTPUTS ONLY ==="
log "Metadata file: $METADATA"
log "Previous versions moved to: $OLD"

# --- sanity: the inputs we are reusing must actually exist -------------------
for f in braycurtis_dm.qza jaccard_dm.qza braycurtis_pcoa.qza jaccard_pcoa.qza \
         merged-table.qza rooted-tree.qza; do
  if [[ ! -s "${MERGED_DIR}/${f}" ]]; then
    log "FATAL: ${f} missing - cannot refresh without it. Aborting, nothing changed."
    exit 1
  fi
done
log "All reusable artefacts present."

# --- move the stale metadata-dependent outputs aside ------------------------
for f in braycurtis_emperor.qzv jaccard_emperor.qzv merged-table-summary.qzv \
         braycurtis_Study_ID_permanova.qzv \
         braycurtis_Intervention_Class_permanova.qzv \
         braycurtis_PCoA_Group_permanova.qzv; do
  [[ -e "${MERGED_DIR}/${f}" ]] && mv "${MERGED_DIR}/${f}" "$OLD/" && log "  archived ${f}"
done

# --- 1. table summary -------------------------------------------------------
log "[1/3] Re-summarising merged table with corrected metadata"
$QIIME feature-table summarize \
  --i-table "${MERGED_DIR}/merged-table.qza" \
  --m-sample-metadata-file "$METADATA" \
  --o-visualization "${MERGED_DIR}/merged-table-summary.qzv"

# --- 2. emperor plots (reusing existing PCoA coordinates) -------------------
# Only the plot is rebuilt; the ordination itself is unchanged because it
# depends on the distance matrix, not the metadata.
for METRIC in braycurtis jaccard; do
  log "[2/3] Re-plotting ${METRIC} PCoA with corrected metadata"
  $QIIME emperor plot \
    --i-pcoa "${MERGED_DIR}/${METRIC}_pcoa.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "${MERGED_DIR}/${METRIC}_emperor.qzv"
done

# --- 3. PERMANOVA -----------------------------------------------------------
# Study_ID is the one that genuinely changes: Ma_2024 is now two studies
# (Ma_2024 = the sequenced FMT trial; the laying-hen trial has no sequencing),
# so the number of batch levels differs from the previous run.
for COL in Study_ID Intervention_Class PCoA_Group; do
  log "[3/3] PERMANOVA on ${COL}"
  $QIIME diversity beta-group-significance \
    --i-distance-matrix "${MERGED_DIR}/braycurtis_dm.qza" \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$COL" \
    --p-method permanova \
    --o-visualization "${MERGED_DIR}/braycurtis_${COL}_permanova.qzv" \
    || log "  (skipped ${COL} - column may be absent or single-level)"
done

log "=== DONE ==="
log "Refreshed, open at https://view.qiime2.org :"
log "  ${MERGED_DIR}/braycurtis_emperor.qzv"
log "  ${MERGED_DIR}/jaccard_emperor.qzv"
log "  ${MERGED_DIR}/braycurtis_Study_ID_permanova.qzv   <- changed: Ma_2024 split"
log "  ${MERGED_DIR}/braycurtis_Intervention_Class_permanova.qzv"
log "Pre-correction versions kept in: ${OLD}"
