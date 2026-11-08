#!/bin/bash
#SBATCH --job-name=pdi_filter
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=/scratch/user/%u/pdi/logs/filter_%j.log
# =============================================================================
# 04_filter_and_export.sh
#
# Applies rare-ASV and contaminant filtering to the merged feature table,
# removes low-depth samples, and exports the final table for downstream
# analysis in R.
#
# Usage (after stage 03):
#   sbatch 04_filter_and_export.sh
#
# Filtering is applied in a fixed order: rare-ASV filter, then contaminant
# filter, then low-depth sample filter. Feature-level filters are applied
# before the sample-depth filter because removing rare and contaminant
# features changes each sample's total read count; filtering samples by
# depth first would evaluate them against a depth that includes reads
# subsequently discarded.
# =============================================================================

# Under sbatch, SLURM copies this script to a spool directory and executes
# it from there, so $(dirname "$0") does not point at the scripts directory.
# SLURM_SUBMIT_DIR is the directory sbatch was invoked from, which is where
# the repository actually lives; this falls back to the script's own
# directory for interactive use on the login node.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
if [[ ! -f "${SCRIPT_DIR}/00_config.sh" ]]; then
  echo "FATAL: 00_config.sh not found in ${SCRIPT_DIR}" >&2
  echo "Submit with:  cd /scratch/user/$USER/pdi/scripts && sbatch $(basename "$0")" >&2
  exit 1
fi
source "${SCRIPT_DIR}/00_config.sh"

if [[ -z "${PDI_ROOT:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh loaded but PDI_ROOT/log() are unset." >&2
  exit 1
fi
activate_qiime

cd "$MERGED_DIR"

# =============================================================================
# STEP 1 - DERIVE THE RARE-ASV THRESHOLD
#
# The filtering rule: the mean per-sample sequencing depth, multiplied by a
# fixed factor (RARE_ASV_MULTIPLIER, set in 00_config.sh) and rounded to the
# nearest integer, becomes the minimum total-frequency threshold for an ASV
# to be retained.
#
# The threshold is computed programmatically from the merged feature-table
# summary rather than read manually from a visualization, so that the value
# used is reproducible and recorded in the log rather than transcribed by
# hand.
# =============================================================================
log "[1/6] Deriving rare-ASV min-frequency threshold"

if [[ ! -s "${MERGED_DIR}/exported-merged-table/feature-table.biom" ]]; then
  $QIIME tools export \
    --input-path "${MERGED_DIR}/merged-table.qza" \
    --output-path "${MERGED_DIR}/exported-merged-table"
fi

biom summarize-table \
  -i "${MERGED_DIR}/exported-merged-table/feature-table.biom" \
  -o "${MERGED_DIR}/merged_table_summary.txt"

head -20 "${MERGED_DIR}/merged_table_summary.txt"

read -r MEAN_DEPTH MEDIAN_DEPTH MIN_FREQ < <(python3 - "$MERGED_DIR/merged_table_summary.txt" "$RARE_ASV_MULTIPLIER" <<'PY'
import re, sys
txt = open(sys.argv[1]).read()
mult = float(sys.argv[2])
def grab(label):
    m = re.search(rf'^\s*{label}:\s*([\d.,]+)\s*$', txt, re.M)
    if not m: sys.exit(f"ERROR: could not parse '{label}' from biom summary")
    return float(m.group(1).replace(',', ''))
mean, median = grab('Mean'), grab('Median')
# Rounds half up to the nearest integer; the threshold is never allowed to
# be zero.
minfreq = max(1, int(mean * mult + 0.5))
print(f"{mean:.2f} {median:.2f} {minfreq}")
PY
)

log "Mean sample depth   : ${MEAN_DEPTH}"
log "Median sample depth : ${MEDIAN_DEPTH}   <- compare these two"
log "min-frequency       : ${MIN_FREQ}   (mean x ${RARE_ASV_MULTIPLIER}, rounded)"
echo "${MIN_FREQ}" > "${MERGED_DIR}/min_frequency_used.txt"

# ---- Skew check -------------------------------------------------------------
# The mean is sensitive to a small number of very high-depth samples in a
# way the median is not. A large gap between the two values indicates the
# threshold above is being driven upward by a minority of high-depth
# samples, which risks discarding genuine low-abundance taxa from
# lower-depth studies. This condition is logged explicitly rather than
# passed through silently.
python3 - "$MEAN_DEPTH" "$MEDIAN_DEPTH" <<'PY'
import sys
mean, med = float(sys.argv[1]), float(sys.argv[2])
if med > 0 and mean / med > 1.5:
    print(f"\n*** SKEW WARNING: mean/median = {mean/med:.2f} ***")
    print("Depth distribution is right-skewed.")
    print("The mean-based threshold may be too aggressive for lower-depth studies.")
    print("Consider a median-based threshold as a sensitivity analysis.\n")
PY

# =============================================================================
# STEP 2 - FILTER RARE ASVs
# --p-min-frequency: an ASV must appear at least this many times, summed
#   across the whole corpus, to be retained.
# --p-min-samples 1: an ASV meeting the frequency threshold in even one
#   sample is kept.
# Singleton and near-singleton ASVs are predominantly residual sequencing
# error. Removing them substantially reduces the feature count while
# stabilizing all downstream steps.
# =============================================================================
if ! already_done "${MERGED_DIR}/table_filt.qza"; then
  log "[2/6] Filtering rare ASVs at min-frequency=${MIN_FREQ}"
  $QIIME feature-table filter-features \
    --i-table "${MERGED_DIR}/merged-table.qza" \
    --p-min-frequency "$MIN_FREQ" \
    --p-min-samples 1 \
    --o-filtered-table "${MERGED_DIR}/table_filt.qza"
fi

# =============================================================================
# STEP 3 - REMOVE CONTAMINANTS AND UNCLASSIFIED ASVs
# --p-exclude mitochondria,chloroplast: host mitochondrial DNA and plant
#   material from ingested feed both amplify with 16S primers. Neither
#   represents gut bacterial signal, and both can be abundant enough to
#   distort relative abundance calculations.
# --p-include p__: retains only ASVs classified to at least phylum level.
#   An ASV that cannot be classified even to phylum cannot contribute to a
#   genus-level index and is excluded.
# =============================================================================
if ! already_done "${MERGED_DIR}/table_filt_contam.qza"; then
  log "[3/6] Removing mitochondria, chloroplast and unclassified ASVs"
  $QIIME taxa filter-table \
    --i-table "${MERGED_DIR}/table_filt.qza" \
    --i-taxonomy "${MERGED_DIR}/taxonomy.qza" \
    --p-include p__ \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table "${MERGED_DIR}/table_filt_contam.qza"

  $QIIME feature-table summarize \
    --i-table "${MERGED_DIR}/table_filt_contam.qza" \
    --m-sample-metadata-file "$METADATA" \
    --o-visualization "${MERGED_DIR}/table_filt_contam_summary.qzv"
fi

# =============================================================================
# STEP 4 - RAREFACTION CURVES (diagnostic; no filtering is applied here)
# Plots observed ASV count against sequencing depth for each sample. A
# curve that flattens indicates the sample has been sequenced deeply enough
# to capture most of its diversity; a curve still rising at the plotted
# maximum indicates undersampling. Used to support the sample-depth cutoff
# applied in step 5.
# =============================================================================
if ! already_done "${MERGED_DIR}/rarefaction_curves.qzv"; then
  log "[4/6] Rarefaction curves"
  MAXD=$(python3 - "${MERGED_DIR}/merged_table_summary.txt" <<'PY'
import re, sys
txt = open(sys.argv[1]).read()
med = float(re.search(r'^\s*Median:\s*([\d.,]+)\s*$', txt, re.M).group(1).replace(',',''))
# The curve is capped at the median depth rather than the corpus maximum.
# Plotting to the true maximum would compress the depth axis for
# lower-depth samples, given a large spread in per-sample sequencing depth
# across the corpus, and make the figure difficult to interpret.
print(int(med))
PY
)
  log "Rarefaction max-depth set to median depth = ${MAXD}"

  # The metric plotted here (observed_features) does not actually use
  # phylogeny - --i-phylogeny is only required by the QIIME2 API signature
  # for this action, not by the metric's math. If 03b_build_phylogeny.sh has
  # not been run yet, rooted-tree.qza will not exist; in that case this step
  # is skipped with a clear note rather than failing, since nothing
  # downstream in this script depends on it either.
  if [[ -s "${MERGED_DIR}/rooted-tree.qza" ]]; then
    $QIIME diversity alpha-rarefaction \
      --i-table "${MERGED_DIR}/table_filt_contam.qza" \
      --i-phylogeny "${MERGED_DIR}/rooted-tree.qza" \
      --m-metadata-file "$METADATA" \
      --p-max-depth "$MAXD" \
      --p-steps 20 \
      --p-metrics observed_features \
      --o-visualization "${MERGED_DIR}/rarefaction_curves.qzv"
  else
    log "  rooted-tree.qza not found (03b_build_phylogeny.sh not run yet) - skipping rarefaction curves."
    log "  This blocks nothing. Run 03b_build_phylogeny.sh then re-run this stage if the curves are needed."
  fi
fi

# =============================================================================
# STEP 5 - DROP LOW-DEPTH SAMPLES, THEN SUBSET REPRESENTATIVE SEQUENCES
#
# SAMPLE_MIN_DEPTH is set to the 25th percentile of the observed per-sample
# read-count distribution, computed after feature-level filtering. This
# value is recorded in the registered Methods and is computed here rather
# than fixed as a constant, so it reflects the corpus as filtered. It
# equalizes sampling effort across retained samples while keeping
# approximately 75% of the total.
# =============================================================================
if ! already_done "${MERGED_DIR}/table_final.qza"; then
  log "[5/6] Dropping low-depth samples"

  $QIIME tools export \
    --input-path "${MERGED_DIR}/table_filt_contam.qza" \
    --output-path "${MERGED_DIR}/exported-filt-contam"
  biom summarize-table \
    -i "${MERGED_DIR}/exported-filt-contam/feature-table.biom" \
    -o "${MERGED_DIR}/filt_contam_summary.txt"

  SAMPLE_MIN_DEPTH=$(python3 - "${MERGED_DIR}/filt_contam_summary.txt" <<'PY'
import re, sys
lines = open(sys.argv[1]).read().splitlines()
# Per-sample counts appear after the 'Counts/sample detail:' header line.
try:
    i = next(k for k,l in enumerate(lines) if 'Counts/sample detail' in l)
except StopIteration:
    sys.exit("ERROR: no per-sample detail in biom summary")
vals = []
for l in lines[i+1:]:
    m = re.search(r':\s*([\d.,]+)\s*$', l)
    if m: vals.append(float(m.group(1).replace(',','')))
vals.sort()
if not vals: sys.exit("ERROR: parsed zero samples")
print(int(vals[int(len(vals)*0.25)]))
PY
)
  log "Sample min depth (25th percentile) = ${SAMPLE_MIN_DEPTH}"
  echo "${SAMPLE_MIN_DEPTH}" > "${MERGED_DIR}/sample_min_depth_used.txt"

  $QIIME feature-table filter-samples \
    --i-table "${MERGED_DIR}/table_filt_contam.qza" \
    --p-min-frequency "$SAMPLE_MIN_DEPTH" \
    --o-filtered-table "${MERGED_DIR}/table_final.qza"

  # Subsets representative sequences to match the final filtered table, so
  # the sequence file and the feature table describe an identical set of
  # ASVs.
  $QIIME feature-table filter-seqs \
    --i-data "${MERGED_DIR}/merged-rep-seqs.qza" \
    --i-table "${MERGED_DIR}/table_final.qza" \
    --o-filtered-data "${MERGED_DIR}/rep_seqs_final.qza"

  $QIIME feature-table summarize \
    --i-table "${MERGED_DIR}/table_final.qza" \
    --m-sample-metadata-file "$METADATA" \
    --o-visualization "${MERGED_DIR}/table_final_summary.qzv"
fi

# =============================================================================
# STEP 6 - EXPORT FOR DOWNSTREAM ANALYSIS IN R
# Downstream steps (batch correction, index construction, model fitting)
# run in R and require the feature table outside QIIME 2's .qza format. The
# BIOM file with taxonomy attached is the primary output; the TSV is a
# human-readable equivalent. Header columns are renamed to match the
# convention expected by biom add-metadata.
# =============================================================================
log "[6/6] Exporting final artefacts"
EXPORT_DIR="${MERGED_DIR}/final_exported"
mkdir -p "$EXPORT_DIR"

$QIIME tools export --input-path "${MERGED_DIR}/table_final.qza"    --output-path "$EXPORT_DIR"
$QIIME tools export --input-path "${MERGED_DIR}/rep_seqs_final.qza" --output-path "$EXPORT_DIR"

cp "${MERGED_DIR}/taxa/taxonomy.tsv" "${EXPORT_DIR}/taxonomy.tsv"
sed -i -e '1 s/Feature ID/#OTUID/' -e '1 s/Taxon/taxonomy/' -e '1 s/Confidence/confidence/' \
  "${EXPORT_DIR}/taxonomy.tsv"

biom add-metadata \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${EXPORT_DIR}/feature-table_w_tax.biom" \
  --observation-metadata-fp "${EXPORT_DIR}/taxonomy.tsv" \
  --observation-header OTUID,taxonomy,confidence \
  --sc-separated taxonomy

biom convert \
  -i "${EXPORT_DIR}/feature-table_w_tax.biom" \
  -o "${EXPORT_DIR}/feature-table_w_tax.tsv" \
  --to-tsv \
  --header-key taxonomy

log "=== STAGE 04 COMPLETE ==="
log "Thresholds used (record these in Methods):"
log "  rare-ASV min-frequency = $(cat "${MERGED_DIR}/min_frequency_used.txt")"
log "  sample min depth       = $(cat "${MERGED_DIR}/sample_min_depth_used.txt")"
