#!/bin/bash
# =============================================================================
# 01g_subsample_wassie.sh
#
# Downsamples Wassie_2021 to 150,000 read pairs prior to denoising.
#
# Usage (login node - I/O bound but no internet access required):
#   bash 01g_subsample_wassie.sh              # dry run - reports the plan
#   bash 01g_subsample_wassie.sh --apply      # performs the subsampling
#
# Rationale: this study was sequenced at approximately 22 million reads per
# sample, roughly 200x the depth of the rest of the corpus (50,000-150,000
# reads per sample). Two consequences follow from this depth disparity:
#
#   1. The rare-ASV filter (stage 04) computes its threshold as mean
#      per-sample depth x 0.001, applied across the merged feature table.
#      Because the mean is outlier-sensitive, one high-depth study raises
#      the threshold applied to every other study.
#
#   2. DADA2's detection sensitivity scales with sequencing depth
#      (approximately 4/depth). At 22 million reads, ASVs are resolved down
#      to a relative abundance the rest of the corpus cannot detect at its
#      typical depth. Since Study_ID is also the batch-correction variable
#      in the downstream MMUPHin step, this creates a confound between
#      sequencing depth and batch identity specific to this study.
#
# Subsampling to 150,000 read pairs brings this study into the same depth
# range as the rest of the corpus.
#
# Subsampling is applied to the raw reads, before QIIME 2 import and before
# denoising - not to the resulting ASV table. Subsampling after denoising
# would not address point 2 above, because DADA2 would already have
# inferred ASVs at a sensitivity no other study in the corpus can match.
#
# The full-depth data is retained in raw_data_full/Wassie_2021/ and is not
# modified by this script, to support a full-depth sensitivity analysis
# alongside the subsampled result.
#
# Read pairing is preserved by using an identical random seed for the
# forward and reverse read files; different seeds would select different
# reads from each file and break mate pairing, causing DADA2's merge step
# to fail.
# =============================================================================

source "$(dirname "$0")/00_config.sh"

if [[ -z "${LOG_DIR:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh did not load correctly." >&2
  exit 1
fi

STUDY="Wassie_2021"
DEPTH=150000          # read PAIRS to keep
SEED=100              # fixed seed - reproducible, and identical for R1/R2
SRC="${PDI_ROOT}/raw_data_full/${STUDY}"
DST="${RAW_DIR}/${STUDY}"

APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1

# ---- seqtk ------------------------------------------------------------------
if ! command -v seqtk >/dev/null 2>&1; then
  echo "seqtk not found. Install it into the pdi env first (LOGIN NODE):"
  echo "  module load Miniconda3/23.10.0-1"
  echo "  source activate /scratch/user/$USER/.conda/envs/pdi"
  echo "  conda install -c bioconda seqtk -y"
  exit 1
fi
log "seqtk: $(command -v seqtk)"

[[ -d "$SRC" ]] || { log "ERROR: full-depth source not found at ${SRC}"; exit 1; }

n_src=$(find "$SRC" -maxdepth 1 -name '*_1.fastq.gz' | wc -l)
log "source runs: ${n_src}  (${SRC})"
log "target      : ${DST}"
log "depth       : ${DEPTH} read pairs   seed: ${SEED}"
[[ "$APPLY" -eq 1 ]] && log "APPLY MODE" || log "DRY RUN - add --apply to execute"
echo

mkdir -p "$DST"
done_n=0; skip_n=0

for r1 in "$SRC"/*_1.fastq.gz; do
  [[ -e "$r1" ]] || continue
  acc=$(basename "$r1" _1.fastq.gz)
  r2="${SRC}/${acc}_2.fastq.gz"
  if [[ ! -e "$r2" ]]; then
    log "  WARN ${acc}: no mate, skipping"
    continue
  fi

  o1="${DST}/${acc}_1.fastq.gz"
  o2="${DST}/${acc}_2.fastq.gz"

  # resumable: a completed pair is left alone
  if [[ -s "$o1" && -s "$o2" ]]; then
    log "  SKIP ${acc} (already subsampled)"
    skip_n=$((skip_n+1)); continue
  fi

  log "  ${acc}"
  if [[ "$APPLY" -eq 1 ]]; then
    # SAME SEED on both mates - this is the part that must not change
    seqtk sample -s"$SEED" "$r1" "$DEPTH" | gzip > "$o1"
    seqtk sample -s"$SEED" "$r2" "$DEPTH" | gzip > "$o2"

    # a subsampled pair with unequal read counts is unusable - fail loudly now
    c1=$(zcat "$o1" | wc -l); c2=$(zcat "$o2" | wc -l)
    if [[ "$c1" -ne "$c2" ]]; then
      log "  ERROR ${acc}: R1 has $((c1/4)) reads, R2 has $((c2/4)) - PAIRING BROKEN"
      rm -f "$o1" "$o2"
      exit 1
    fi
    log "     -> $((c1/4)) read pairs"
  fi
  done_n=$((done_n+1))
done

echo
if [[ "$APPLY" -eq 0 ]]; then
  log "DRY RUN complete: ${done_n} runs would be subsampled, ${skip_n} already done."
  log "Re-run with --apply."
else
  log "Subsampled ${done_n} runs (${skip_n} already present)."
  log "Full-depth copy untouched at ${SRC}"
  echo
  log "Record for Methods:"
  log "  Wassie_2021 subsampled to ${DEPTH} read pairs/sample, seqtk seed ${SEED},"
  log "  applied to raw reads prior to denoising. Full-depth data retained for"
  log "  sensitivity analysis."
fi
