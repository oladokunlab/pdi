#!/bin/bash
#SBATCH --job-name=pdi_ordination
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/user/%u/pdi/qc_reports/ordination_%j.log
#SBATCH --error=/scratch/user/%u/pdi/qc_reports/ordination_%j.err
# Add your allocation if sbatch refuses the job:
##SBATCH --account=<your_account>
# =============================================================================
# Rebuilds the class ordination and PERMANOVA on the corrected eight-class
# labels, replacing the artefacts built on Intervention_Class and PCoA_Group.
#
# Run either way:
#   sbatch run_ordination.sh          (preferred)
#   bash   run_ordination.sh          (interactive node; not a login node,
#                                      the pairwise PERMANOVA is heavy)
#
# Idempotent: every step is skipped if its output already exists, so a
# re-run after a timeout picks up where it stopped. Delete an output to
# force that step to run again.
# =============================================================================
set -euo pipefail

ROOT=/scratch/user/$USER/pdi
MERGED=$ROOT/merged
META=$ROOT/scripts/GDI_metadata_qiime2_class.tsv
RSCRIPT=/sw/eb/sw/R/4.5.2-gfbf-2025b/bin/Rscript
export R_LIBS_USER=/scratch/user/$USER/R/library
export PDI_ROOT=$ROOT
export CLDN1_MODE=merge

say () { echo "[$(date '+%F %T')] $*"; }
need () { [ -s "$1" ] || { echo "MISSING: $1" >&2; exit 1; }; }

say "=== PDI ordination rebuild ==="
need "$META"
need "$MERGED/braycurtis_dm.qza"
need "$MERGED/jaccard_dm.qza"
need "$MERGED/braycurtis_pcoa.qza"
need "$MERGED/jaccard_pcoa.qza"
say "inputs present"

# ---------------------------------------------------------------------------
# 1. QIIME: PERMANOVA and Emperor on the corrected labels
#
# These run on the merged table of all 3,037 runs across 50 studies, which
# mixes V3-V4 with V4, V3-V5 and full-length amplicons, and they permute
# samples freely with no blocking within study. Both facts inflate
# significance and must be stated wherever the numbers are used. The
# defensible test is step 2.
# ---------------------------------------------------------------------------
module purge
module load QIIME2/2024.10-Amplicon 2>/dev/null || module load QIIME2
say "QIIME: $(qiime --version 2>&1 | head -1)"

cd "$MERGED"
for M in braycurtis jaccard; do
  OUT=${M}_Class_permanova.qzv
  if [ -s "$OUT" ]; then
    say "skip $OUT (exists)"
  else
    say "PERMANOVA on $M, 9999 permutations, pairwise"
    qiime diversity beta-group-significance \
      --i-distance-matrix ${M}_dm.qza \
      --m-metadata-file "$META" \
      --m-metadata-column Class \
      --p-method permanova \
      --p-permutations 9999 \
      --p-pairwise \
      --o-visualization "$OUT"
    say "wrote $OUT"
  fi

  EMP=${M}_Class_emperor.qzv
  if [ -s "$EMP" ]; then
    say "skip $EMP (exists)"
  else
    qiime emperor plot \
      --i-pcoa ${M}_pcoa.qza \
      --m-metadata-file "$META" \
      --o-visualization "$EMP"
    say "wrote $EMP"
  fi
done

# Export the PERMANOVA tables to plain text so the numbers are readable
# without opening a qzv, and so they can go straight into the deposit.
mkdir -p "$ROOT/r_analysis/verification/qiime_class_permanova"
for M in braycurtis jaccard; do
  D="$ROOT/r_analysis/verification/qiime_class_permanova/$M"
  rm -rf "$D"
  qiime tools export --input-path ${M}_Class_permanova.qzv --output-path "$D"
  say "exported $M -> $D"
done

# Retire the superseded artefacts rather than deleting them, so the old
# numbers stay traceable to the vocabulary that produced them.
mkdir -p "$MERGED/superseded_vocab"
for F in braycurtis_Intervention_Class_permanova.qzv \
         braycurtis_PCoA_Group_permanova.qzv; do
  [ -e "$MERGED/$F" ] && mv "$MERGED/$F" "$MERGED/superseded_vocab/" && say "retired $F"
done
cat > "$MERGED/superseded_vocab/README.txt" <<'EOF'
Retired 2026-10-23. These were computed on Intervention_Class and PCoA_Group.

Intervention_Class has six levels that conflate insults with additives: its
Dietary_Intervention level holds 466 additive runs plus 30 toxin challenge and
22 dietary challenge, and its Infectious_Challenge level holds 94 infectious,
59 toxin, 50 challenged-and-treated and 15 additive.

PCoA_Group buckets seven levels as "Other <something>", which the
classification specification forbids.

Replaced by braycurtis_Class_permanova.qzv and jaccard_Class_permanova.qzv,
computed on the eight-class scheme in scripts/GDI_metadata_qiime2_class.tsv.
EOF

# ---------------------------------------------------------------------------
# 2. R: the defensible test - Aitchison distance on the V3-V4 batch-adjusted
#    table, permutations blocked within study, betadisper alongside.
# ---------------------------------------------------------------------------
module purge
if [ -s "$ROOT/r_analysis/gdi/pcoa_strata.tsv" ]; then
  say "skip 55_class_ordination.R (pcoa_strata.tsv exists)"
else
  say "R: 55_class_ordination.R"
  "$RSCRIPT" "$ROOT/scripts/R_scripts/55_class_ordination.R"
fi

# ---------------------------------------------------------------------------
# 3. Python: rebuild every figure, including the new ordination panel
# ---------------------------------------------------------------------------
module purge
module load Anaconda3/2025.12-2
say "Python: make_figures.py"
cd "$ROOT"
python make_figures.py \
  --root   "$ROOT" \
  --labels "$ROOT/intervention_classification_FINAL.csv" \
  --outdir "$ROOT/figures"

say "=== DONE ==="
echo
echo "---- blocked-permutation PERMANOVA (the one to report) ----"
sed -n '1,60p' "$ROOT/r_analysis/verification/55_class_ordination.txt" || true
echo
echo "---- QIIME class PERMANOVA, free permutation, region-mixed ----"
for M in braycurtis jaccard; do
  echo "[$M]"
  head -20 "$ROOT/r_analysis/verification/qiime_class_permanova/$M/index.html" 2>/dev/null \
    | sed 's/<[^>]*>//g' | grep -v '^\s*$' | head -10 || true
done
