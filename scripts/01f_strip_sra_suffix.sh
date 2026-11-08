#!/bin/bash
# =============================================================================
# 01f_strip_sra_suffix.sh - remove the ".sra" infix from every FASTQ filename
#
# HOW TO RUN: login node. DRY RUN BY DEFAULT.
#   bash 01f_strip_sra_suffix.sh            # show what would change
#   bash 01f_strip_sra_suffix.sh --apply    # actually rename
#
# WHY: fastq-dump/fasterq-dump leave the ".sra" container name inside the
# output filename:
#     SRR26348459.sra_1.fastq.gz   instead of   SRR26348459_1.fastq.gz
#     SRR12178618.sra.fastq.gz     instead of   SRR12178618.fastq.gz
#
# The manifest builder (01_build_manifests.sh) derives the QIIME 2 sample-id
# from the filename. With the ".sra" infix left in place, the sample-id
# becomes "SRR26348459.sra", which does not match the plain accession in the
# metadata file, causing the feature-table merge in stage 03 to fail on
# unmatched samples. 01_build_manifests.sh strips this infix defensively, but
# normalizing filenames on disk keeps naming consistent across the corpus.
#
# Handles all three shapes:
#   *.sra_1.fastq.gz -> *_1.fastq.gz
#   *.sra_2.fastq.gz -> *_2.fastq.gz
#   *.sra.fastq.gz   -> *.fastq.gz
# Also covers uncompressed .fastq in case anything slipped through 01a.
#
# NEVER OVERWRITES: if the target name already exists, it reports a collision
# and skips, so a half-renamed folder cannot silently lose a file.
# =============================================================================

source "$(dirname "$0")/00_config.sh"

if [[ -z "${LOG_DIR:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh did not load correctly." >&2
  exit 1
fi

APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1
[[ "$APPLY" -eq 1 ]] && log "APPLY MODE - files will be renamed" \
                     || log "DRY RUN - nothing will change. Add --apply to execute."
echo

total=0; collided=0; studies=0

for d in "${RAW_DIR}"/*/ "${PDI_ROOT}"/raw_data_full/*/; do
  [[ -d "$d" ]] || continue
  study=$(basename "${d%/}")
  n=0

  # -maxdepth 1: only the study folder itself, not any nested leftovers
  while IFS= read -r f; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f"); dir=$(dirname "$f")
    new="${base/.sra_/_}"      # .sra_1.fastq.gz -> _1.fastq.gz
    new="${new/.sra./.}"       # .sra.fastq.gz   -> .fastq.gz
    [[ "$new" == "$base" ]] && continue

    if [[ -e "${dir}/${new}" ]]; then
      log "  COLLISION ${study}: ${base} -> ${new} already exists, SKIPPING"
      collided=$((collided+1)); continue
    fi

    [[ "$APPLY" -eq 1 ]] && mv -n "$f" "${dir}/${new}"
    n=$((n+1))
  done < <(find "$d" -maxdepth 1 -name '*.sra*' \( -name '*.fastq' -o -name '*.fastq.gz' \) 2>/dev/null)

  if [[ "$n" -gt 0 ]]; then
    log "${study}: ${n} files"
    total=$((total+n)); studies=$((studies+1))
  fi
done

echo
log "${total} files across ${studies} studies"
[[ "$collided" -gt 0 ]] && log "${collided} collisions skipped - resolve these by hand"
if [[ "$APPLY" -eq 0 ]]; then
  log "Re-run with --apply to rename."
else
  log "Done. Verify with:  find ${RAW_DIR} -name '*.sra*' | wc -l   (expect 0)"
fi
