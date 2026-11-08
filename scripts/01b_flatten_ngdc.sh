#!/bin/bash
# =============================================================================
# 01b_flatten_ngdc.sh - flatten nested NGDC/NMDC downloads into study folders
#
# HOW TO RUN: login node. Dry-run FIRST, then apply.
#   bash 01b_flatten_ngdc.sh              # dry run, changes nothing
#   bash 01b_flatten_ngdc.sh --apply      # actually move/rename
#
# SRA-toolkit downloads write files flat into the study folder. NGDC/NMDC
# downloads retrieved via `wget -r`,
# which mirrors the server's directory tree instead:
#
#   Du_2026/CRA024363/CRR1716066/CRR1716066_r1.fq.gz
#   ^study  ^project  ^run        ^file
#
# Two things break because of that:
#   1. DEPTH. Both 01a_normalize and 01_build_manifests use -maxdepth 1 /
#      non-recursive globs, so they never see files three levels down. The
#      files aren't rejected - they're simply invisible, which is the worst
#      kind of failure.
#   2. NAMING. NGDC uses lowercase _r1/_r2 and the short .fq.gz extension.
#      Stage 01's globs are case-sensitive and match _R1/_1 with .fastq.gz.
#      So even flattened, they'd still be missed.
#
# This script fixes both: moves every FASTQ up into <study>/, and renames
#   <run>_r1.fq.gz -> <run>_1.fastq.gz
#   <run>_r2.fq.gz -> <run>_2.fastq.gz
# so the result is indistinguishable from an SRA-toolkit download.
#
# SAFE BY DEFAULT: dry-run unless --apply is passed. Never overwrites an
# existing target - it reports a collision and skips instead.
# =============================================================================

source "$(dirname "$0")/00_config.sh"

if [[ -z "${LOG_DIR:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh did not load correctly." >&2
  exit 1
fi

APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1

if [[ "$APPLY" -eq 0 ]]; then
  log "DRY RUN - nothing will be changed. Re-run with --apply to execute."
else
  log "APPLY MODE - files will be moved and renamed."
fi
echo

moved=0; collided=0; studies_touched=0

tail -n +2 "$STUDY_TABLE" | while IFS=$'\t' read -r study n_runs region platform layout pf pr status; do
  [[ -z "$study" ]] && continue
  sdir="${RAW_DIR}/${study}"
  [[ -d "$sdir" ]] || continue

  # Only act on studies that actually have FASTQ files nested below depth 1.
  nested=$(find "$sdir" -mindepth 2 \( -iname '*.fq' -o -iname '*.fq.gz' \
             -o -iname '*.fastq' -o -iname '*.fastq.gz' \) 2>/dev/null | wc -l)
  [[ "$nested" -eq 0 ]] && continue

  log "=== ${study}: ${nested} nested FASTQ found (expect ${n_runs} runs) ==="
  studies_touched=$((studies_touched+1))

  find "$sdir" -mindepth 2 \( -iname '*.fq' -o -iname '*.fq.gz' \
       -o -iname '*.fastq' -o -iname '*.fastq.gz' \) | sort | while read -r f; do
    base=$(basename "$f")

    # Normalize extension: .fq -> .fastq, .fq.gz -> .fastq.gz
    new="${base/.fq.gz/.fastq.gz}"
    new="${new/.fq/.fastq}"
    # Normalize mate marker: lowercase _r1/_r2 -> _1/_2 (what stage 01 matches)
    new="${new/_r1./_1.}"
    new="${new/_r2./_2.}"
    new="${new/_R1./_1.}"
    new="${new/_R2./_2.}"

    target="${sdir}/${new}"

    if [[ -e "$target" ]]; then
      log "  COLLISION: ${base} -> ${new} already exists at top level, SKIPPING"
      collided=$((collided+1))
      continue
    fi

    if [[ "$APPLY" -eq 1 ]]; then
      mv -n "$f" "$target"
    fi
    log "  ${base}  ->  ${new}"
    moved=$((moved+1))
  done

  if [[ "$APPLY" -eq 1 ]]; then
    # Remove the now-empty mirrored directory tree. -type d -empty only
    # deletes directories with nothing left in them, so this cannot take
    # data with it.
    find "$sdir" -mindepth 1 -type d -empty -delete 2>/dev/null || true
    remaining=$(find "$sdir" -maxdepth 1 -iname '*.fastq.gz' | wc -l)
    log "  ${study}: now ${remaining} files at top level (expect $((n_runs * 2)) if paired, ${n_runs} if single)"
  fi
  echo
done

echo
if [[ "$APPLY" -eq 0 ]]; then
  log "DRY RUN complete. Nothing changed."
  log "If the renames above look right, re-run:  bash $(basename "$0") --apply"
else
  log "APPLY complete."
  log "Next: re-run 01a_normalize_raw_data.sh (catches anything left), then 01_build_manifests.sh"
fi
