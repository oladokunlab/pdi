#!/bin/bash
#SBATCH --job-name=pdi_picrust2
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --partition=medium
#SBATCH --output=/scratch/user/%u/pdi/logs/picrust2_%j.log
# =============================================================================
# 05_picrust2.sh
#
# Predicts functional gene content (KEGG orthologs and MetaCyc pathways) from
# the 16S ASVs, so the community can be analysed by what it can DO rather than
# by which genera are present.
#
# SUBMIT:  cd /scratch/user/$USER/pdi/scripts && sbatch 05_picrust2.sh
# CHECK:   squeue -u $USER
#
# -----------------------------------------------------------------------------
# WHY THIS IS WORTH RUNNING GIVEN EVERYTHING THAT CAME BEFORE
# -----------------------------------------------------------------------------
# The run-level work established that classifiers were memorising study
# identity rather than learning a dysbiosis signature. Part of the reason is
# technical: different primers, regions and pipelines recover different
# GENERA, so genus composition carries a study fingerprint.
#
# Function may not. Two studies can detect entirely different butyrate
# producers - different names, same job - and at the pathway level those
# converge. If the study-specific component of the signal is largely
# taxonomic, a functional representation should transfer better.
#
# That is a testable hypothesis, not an expectation. If pathways behave the
# same as genera, the negative result is worth reporting: it closes the
# obvious reviewer question of whether a functional representation would have
# worked.
#
# -----------------------------------------------------------------------------
# WHAT PICRUSt2 IS AND IS NOT
# -----------------------------------------------------------------------------
# It places each ASV on a reference phylogeny of genomes with known gene
# content, then infers what genes that ASV's organism likely carries. The
# output is a PREDICTION, not a measurement. Accuracy is weaker for
# environments under-represented in the reference set, and chicken gut is
# less well covered than human gut or soil.
#
# The NSTI value quantifies this per ASV: the phylogenetic distance to the
# nearest reference genome. Higher means the prediction rests on a more
# distant relative. ASVs above NSTI 2.0 are excluded by default; the run
# summary reports how many were lost, which is the honest measure of how far
# these predictions can be trusted here.
#
# -----------------------------------------------------------------------------
# RESOURCES
# -----------------------------------------------------------------------------
# Roughly 35,000 ASVs. Sequence placement (EPA-ng) is the expensive step in
# both time and memory. 24 cores and 180 GB with a 12 hour limit is generous;
# if it finishes far short of that, the limits can be lowered for any rerun.
# =============================================================================

set -uo pipefail

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

PIC_DIR="${PDI_ROOT}/picrust2"
mkdir -p "$PIC_DIR" "${PDI_ROOT}/logs"

log "=== PICRUSt2 functional prediction ==="

# --- 1. locate the ASV sequences and feature table ---------------------------
# PICRUSt2 predicts per ASV, so it needs the PRE-agglomeration artefacts - the
# genus table cannot be used. 04_filter_and_export.sh already produced these
# in merged/final_exported, and they remain valid: they derive from sequences
# and counts, neither of which the metadata corrections touched (those changed
# study names, arm names and ages, while the feature table is keyed by run
# accession).
EXPORT_DIR="${MERGED_DIR}/final_exported"
log "[1/4] Locating ASV sequences and feature table"
if [[ -s "${EXPORT_DIR}/dna-sequences.fasta" && -s "${EXPORT_DIR}/feature-table.biom" ]]; then
  log "  using existing exports in ${EXPORT_DIR}"
else
  log "  exports not found - regenerating from the QIIME2 artefacts"
  EXPORT_DIR="${PIC_DIR}/input"
  mkdir -p "$EXPORT_DIR"
  activate_qiime
  $QIIME tools export --input-path "${MERGED_DIR}/rep_seqs_final.qza" --output-path "$EXPORT_DIR"
  $QIIME tools export --input-path "${MERGED_DIR}/table_final.qza"    --output-path "$EXPORT_DIR"
  log "  exported to ${EXPORT_DIR}"
fi

for f in dna-sequences.fasta feature-table.biom; do
  if [[ ! -s "${EXPORT_DIR}/${f}" ]]; then
    log "FATAL: ${EXPORT_DIR}/${f} missing or empty after export. Aborting."
    exit 1
  fi
done
N_ASV=$(grep -c '^>' "${EXPORT_DIR}/dna-sequences.fasta")
log "  ${N_ASV} ASV sequences ready"

# --- 2. load PICRUSt2 --------------------------------------------------------
log "[2/4] Loading PICRUSt2"
module purge
module load GCC/12.3.0 OpenMPI/4.1.5 PICRUSt2/2.6.2
if ! command -v picrust2_pipeline.py >/dev/null 2>&1; then
  log "FATAL: picrust2_pipeline.py not on PATH after module load. Aborting."
  exit 1
fi
log "  $(picrust2_pipeline.py --version 2>&1 | head -1)"

# --- 3. run the pipeline -----------------------------------------------------
# --stratified is deliberately NOT used: it produces per-taxon contributions
# to each pathway, which multiplies output size enormously and is not needed
# for the community-level analysis planned here.
OUT="${PIC_DIR}/out"
if [[ -d "$OUT" ]]; then
  log "  removing previous output directory so the run starts clean"
  rm -rf "$OUT"
fi

log "[3/4] Running picrust2_pipeline.py (this is the long step)"
picrust2_pipeline.py \
  -s "${EXPORT_DIR}/dna-sequences.fasta" \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "$OUT" \
  -p "${SLURM_CPUS_PER_TASK:-24}" \
  --in_traits EC,KO \
  --verbose

# --- 4. decompress the tables R will read, and report NSTI ------------------
log "[4/4] Decompressing output tables"
for f in "${OUT}"/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
         "${OUT}"/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
         "${OUT}"/pathways_out/path_abun_unstrat.tsv.gz \
         "${OUT}"/*/weighted_nsti.tsv.gz; do
  [[ -e "$f" ]] && gunzip -kf "$f" && log "  $(basename "$f")"
done

NSTI="${OUT}/KO_metagenome_out/weighted_nsti.tsv"
if [[ -s "$NSTI" ]]; then
  log ""
  log "NSTI - how far each sample's ASVs sit from a reference genome:"
  awk 'NR>1 {s+=$2; n++; if($2>m) m=$2} END {
    printf "  mean %.3f, max %.3f, across %d samples\n", s/n, m, n }' "$NSTI"
  log "  Below ~0.15 is usually considered reliable for well-characterised"
  log "  environments. Chicken gut is under-represented in the reference set,"
  log "  so a higher value here is expected and must be reported as a caveat."
fi

log ""
log "=== DONE ==="
log "Tables for R:"
log "  ${OUT}/KO_metagenome_out/pred_metagenome_unstrat.tsv     (KEGG orthologs)"
log "  ${OUT}/pathways_out/path_abun_unstrat.tsv                (MetaCyc pathways)"
log "  ${NSTI}                                                   (prediction quality)"
