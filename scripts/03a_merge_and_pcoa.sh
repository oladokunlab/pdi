#!/bin/bash
#SBATCH --job-name=pdi_pcoa
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=/scratch/user/%u/pdi/logs/pcoa_%j.log
# =============================================================================
# 03a_merge_and_pcoa.sh - merge the corpus and produce a PCoA, fast
#
# HOW TO RUN:  cd /scratch/user/$USER/pdi/scripts && sbatch 03a_merge_and_pcoa.sh
#
# WHY THIS EXISTS SEPARATELY FROM 03_merge_and_taxonomy.sh:
# The full stage 03 is dominated by classify-sklearn against the full-length
# SILVA classifier, which at corpus scale takes 4-12 hours on its own. This does
# not need taxonomy to look at community structure - PCoA needs only the merged
# feature table. So this script does merge -> Bray-Curtis -> PCoA and stops,
# giving me a figure in about an hour instead of overnight.
#
# Bray-Curtis rather than UniFrac is deliberate here: UniFrac needs a phylogeny,
# and MAFFT + FastTree on this many ASVs is another 2-5 hours. Bray-Curtis needs
# no tree. It is also the more appropriate choice for this corpus - aligning V3-V4, V4
# and PacBio full-length sequences in one alignment is imperfect, so a
# phylogeny-free metric avoids this caveat.
#
# NOTHING HERE IS WASTED: merged-table.qza is the same artefact stage 03 needs,
# so running the full taxonomy pass afterwards reuses it and skips the merge.
#
# WHAT TO LOOK FOR:
#   Coloured by Study_ID  -> if samples cluster by study, that is the technical
#                            batch effect MMUPHin exists to correct. Seeing it
#                            here is EXPECTED, not a failure.
#   Coloured by PCoA_Group -> whether treatment structure is visible at all
#                            underneath the batch effect.
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

# ---------------------------------------------------------------------------
# STEP 1 - MERGE FEATURE TABLES
# 51 studies denoised independently become one table. Valid because ASVs are
# exact sequences: the feature ID is an md5 hash of the sequence, so identical
# sequences from different studies collapse to the same feature automatically.
# error_on_overlapping_sample is intentional - if two studies somehow share a
# sample ID, this produces a hard failure rather than a silent overwrite.
# ---------------------------------------------------------------------------
if ! already_done "${MERGED_DIR}/merged-table.qza"; then
  log "[1/5] Merging per-study tables"
  TABLES=(); SEQS=()
  while IFS=$'\t' read -r study rest; do
    [[ -z "$study" ]] && continue
    t="${WORK_DIR}/${study}/table.qza"; s="${WORK_DIR}/${study}/rep-seqs.qza"
    if [[ -s "$t" && -s "$s" ]]; then
      TABLES+=(--i-tables "$t"); SEQS+=(--i-data "$s")
    else
      log "  WARN ${study}: no denoised output, EXCLUDED from merge"
    fi
  done < <(tail -n +2 "$STUDY_TABLE" | tr -d '\r')
  log "  merging $(( ${#TABLES[@]} / 2 )) studies"

  $QIIME feature-table merge "${TABLES[@]}" \
    --p-overlap-method error_on_overlapping_sample \
    --o-merged-table "${MERGED_DIR}/merged-table.qza"

  $QIIME feature-table merge-seqs "${SEQS[@]}" \
    --o-merged-data "${MERGED_DIR}/merged-rep-seqs.qza"
fi

# ---------------------------------------------------------------------------
# STEP 2 - SUMMARISE
# This .qzv is where stage 04 reads MEAN SAMPLE DEPTH from for the rare-ASV
# threshold (mean x 0.001) in stage 04. It also provides an initial view of
# the whole corpus.
# ---------------------------------------------------------------------------
if ! already_done "${MERGED_DIR}/merged-table-summary.qzv"; then
  log "[2/5] Summarising merged table"
  $QIIME feature-table summarize \
    --i-table "${MERGED_DIR}/merged-table.qza" \
    --m-sample-metadata-file "$METADATA" \
    --o-visualization "${MERGED_DIR}/merged-table-summary.qzv"
fi

# ---------------------------------------------------------------------------
# STEP 3 - RAREFY
# PCoA needs even sampling effort, or samples separate by depth rather than by
# composition. The 25th percentile of the observed per-sample depth is used,
# matching the value specified in the registered Methods, retaining
# approximately 75% of samples.
# ---------------------------------------------------------------------------
if ! already_done "${MERGED_DIR}/rarefied-table.qza"; then
  log "[3/5] Deriving rarefaction depth"
  $QIIME tools export --input-path "${MERGED_DIR}/merged-table.qza" \
                      --output-path "${MERGED_DIR}/exp-merged" >/dev/null
  biom summarize-table -i "${MERGED_DIR}/exp-merged/feature-table.biom" \
                       -o "${MERGED_DIR}/merged_summary.txt"

  DEPTH=$(python3 - "${MERGED_DIR}/merged_summary.txt" <<'PY'
import re, sys
lines = open(sys.argv[1]).read().splitlines()
i = next(k for k,l in enumerate(lines) if 'Counts/sample detail' in l)
v = []
for l in lines[i+1:]:
    m = re.search(r':\s*([\d.,]+)\s*$', l)
    if m: v.append(float(m.group(1).replace(',','')))
v.sort()
print(int(v[int(len(v)*0.25)]))
PY
)
  log "  rarefying to ${DEPTH} (25th percentile)"
  echo "$DEPTH" > "${MERGED_DIR}/pcoa_depth_used.txt"

  $QIIME feature-table rarefy \
    --i-table "${MERGED_DIR}/merged-table.qza" \
    --p-sampling-depth "$DEPTH" \
    --o-rarefied-table "${MERGED_DIR}/rarefied-table.qza"
fi

# ---------------------------------------------------------------------------
# STEP 4 - BRAY-CURTIS + PCoA
# Bray-Curtis measures compositional dissimilarity from abundances alone, no
# phylogeny required. Jaccard added as a presence/absence cross-check: if the
# two disagree, the difference between samples is in WHICH taxa are present
# rather than in how abundant they are.
# ---------------------------------------------------------------------------
for METRIC in braycurtis jaccard; do
  if ! already_done "${MERGED_DIR}/${METRIC}_pcoa.qza"; then
    log "[4/5] ${METRIC}"
    $QIIME diversity beta \
      --i-table "${MERGED_DIR}/rarefied-table.qza" \
      --p-metric "$METRIC" \
      --p-n-jobs "$NCORES" \
      --o-distance-matrix "${MERGED_DIR}/${METRIC}_dm.qza"
    $QIIME diversity pcoa \
      --i-distance-matrix "${MERGED_DIR}/${METRIC}_dm.qza" \
      --o-pcoa "${MERGED_DIR}/${METRIC}_pcoa.qza"
    $QIIME emperor plot \
      --i-pcoa "${MERGED_DIR}/${METRIC}_pcoa.qza" \
      --m-metadata-file "$METADATA" \
      --o-visualization "${MERGED_DIR}/${METRIC}_emperor.qzv"
  fi
done

# ---------------------------------------------------------------------------
# STEP 5 - PERMANOVA
# Quantifies what the PCoA shows. Study_ID tests the batch effect;
# Intervention_Class and PCoA_Group test whether treatment structure survives
# it. Comparing the two R2 values tells me how much of the variance is
# technical versus biological BEFORE any batch correction.
#
# No --p-pairwise on Study_ID: 51 levels makes pairwise testing very slow and
# only the overall statistic is needed here.
# ---------------------------------------------------------------------------
for COL in Study_ID Intervention_Class PCoA_Group; do
  out="${MERGED_DIR}/braycurtis_${COL}_permanova.qzv"
  already_done "$out" && continue
  log "[5/5] PERMANOVA on ${COL}"
  $QIIME diversity beta-group-significance \
    --i-distance-matrix "${MERGED_DIR}/braycurtis_dm.qza" \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$COL" \
    --p-method permanova \
    --o-visualization "$out" || log "  (skipped ${COL})"
done

log "=== DONE ==="
log "Download and open at https://view.qiime2.org :"
log "  ${MERGED_DIR}/braycurtis_emperor.qzv     <- colour by Study_ID, then PCoA_Group"
log "  ${MERGED_DIR}/jaccard_emperor.qzv"
log "  ${MERGED_DIR}/merged-table-summary.qzv"
log "  ${MERGED_DIR}/braycurtis_*_permanova.qzv"
