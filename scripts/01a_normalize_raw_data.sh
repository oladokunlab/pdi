#!/bin/bash
# =============================================================================
# 01a_normalize_raw_data.sh - canonicalize FASTQ naming BEFORE building manifests
#
# HOW TO RUN: login node. Local filesystem operations only, no internet, so
# this is safe as nohup bash - but for very large studies (Wassie_2021 at
# ~22M reads/sample) the gzip step is CPU-bound and can run a while. If Grace's
# login-node policy limits sustained CPU, run this per-study or during an
# interactive salloc session instead of one long nohup for everything.
#
#   nohup bash 01a_normalize_raw_data.sh > logs/01a_normalize.log 2>&1 &
#
# WHY THIS EXISTS: 01_build_manifests.sh only recognises these patterns:
#   *_1.fastq.gz  *_1.fastq  *_R1*.fastq.gz  *_R1*.fastq   (and _2/_R2 mates)
# Nothing else. Files coming from four different sources (SRA-toolkit, ENA,
# NGDC, NMDC) do not all use that naming; .fq.gz, .fq, and uncompressed
# .fastq extensions occur, sometimes mixed within a single study folder.
# Anything not matching the pattern above is INVISIBLE to the manifest
# builder - it won't error, it will just silently not appear in the manifest,
# which is worse than a loud failure. This script fixes that up front.
#
# WHAT THIS DOES, per study folder in raw_data/:
#   1. Renames *.fq.gz          -> *.fastq.gz   (metadata-only, instant)
#   2. Compresses *.fastq       -> *.fastq.gz   (gzip, keeps data)
#   3. Renames+compresses *.fq  -> *.fastq.gz
#   4. Detects and WARNS about (does not silently fix):
#      - name collisions where a normalized target already exists
#      - "bare" accession files with no _1/_2/R1/R2 mate marker in a study
#        whose layout is 'paired' - these are genuinely ambiguous (could be
#        single-end despite the metadata, could be an interleaved paired file,
#        could be a joined/merged-read file); these are flagged for manual review
# rather than resolved automatically.
#
# RESUMABLE: already-correct .fastq.gz files are left untouched; re-running
# this after fixing a flagged study only touches what's still wrong.
# =============================================================================

source "$(dirname "$0")/00_config.sh"

# Belt-and-suspenders: if 00_config.sh failed to load for any reason (wrong
# path, typo, missing file), LOG_DIR and the log() function would be unset and
# everything below would silently do the wrong thing instead of stopping. Fail
# loudly here instead.
if [[ -z "${LOG_DIR:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh did not load correctly." >&2
  echo "Expected it at: $(dirname "$0")/00_config.sh" >&2
  echo "Check that file exists there and re-run." >&2
  exit 1
fi

REPORT="${LOG_DIR}/01a_normalize_report.tsv"
: > "$REPORT"
printf 'study\tissue\tdetail\n' >> "$REPORT"

# Use pigz (parallel gzip) if available - meaningfully faster on the large
# studies (Biagi_2020, Miska_2024, Wassie_2021). Falls back to gzip otherwise.
if command -v pigz >/dev/null 2>&1; then
  GZIP_CMD="pigz -p ${NCORES}"
  log "Using pigz (${NCORES} threads) for compression"
else
  GZIP_CMD="gzip"
  log "pigz not found, using single-threaded gzip (slower on large studies)"
fi

n_renamed=0
n_compressed=0
n_collisions=0
n_ambiguous=0

tail -n +2 "$STUDY_TABLE" | while IFS=$'\t' read -r study n_runs region platform layout pf pr status; do
  [[ -z "$study" ]] && continue
  sdir="${RAW_DIR}/${study}"
  [[ -d "$sdir" ]] || continue

  log "=== ${study} (${layout}) ==="

  # ---- 1. Rename *.fq.gz -> *.fastq.gz (no recompression needed) -----------
  find "$sdir" -maxdepth 1 -iname '*.fq.gz' | while read -r f; do
    target="${f%.fq.gz}.fastq.gz"
    if [[ -e "$target" && "$target" != "$f" ]]; then
      echo -e "${study}\tCOLLISION\t${f} -> ${target} (target already exists)" >> "$REPORT"
      log "  COLLISION: $(basename "$f") -> $(basename "$target") already exists, NOT overwriting"
      continue
    fi
    mv -n "$f" "$target"
    log "  renamed: $(basename "$f") -> $(basename "$target")"
  done

  # ---- 2. Compress *.fastq -> *.fastq.gz ------------------------------------
  find "$sdir" -maxdepth 1 -iname '*.fastq' ! -iname '*.fastq.gz' | while read -r f; do
    target="${f}.gz"
    if [[ -e "$target" ]]; then
      echo -e "${study}\tCOLLISION\t${f} (both .fastq and .fastq.gz exist)" >> "$REPORT"
      log "  COLLISION: $(basename "$f") - both compressed and uncompressed exist, NOT touching"
      continue
    fi
    log "  compressing: $(basename "$f")"
    $GZIP_CMD "$f"
  done

  # ---- 3. Rename+compress *.fq -> *.fastq.gz --------------------------------
  find "$sdir" -maxdepth 1 -iname '*.fq' ! -iname '*.fq.gz' | while read -r f; do
    renamed="${f%.fq}.fastq"
    target="${renamed}.gz"
    if [[ -e "$target" || -e "$renamed" ]]; then
      echo -e "${study}\tCOLLISION\t${f} -> ${target} (target already exists)" >> "$REPORT"
      log "  COLLISION: $(basename "$f") -> target already exists, NOT touching"
      continue
    fi
    mv -n "$f" "$renamed"
    log "  compressing: $(basename "$renamed")"
    $GZIP_CMD "$renamed"
  done

  # ---- 4. Flag ambiguous bare files in paired-layout studies ----------------
  # A "bare" file has no _1/_2/R1/R2 marker at all. In a paired study this is
  # genuinely ambiguous and is flagged for manual review rather than resolved
  # automatically.
  if [[ "$layout" == "paired" ]]; then
    find "$sdir" -maxdepth 1 -iname '*.fastq.gz' \
      ! -iname '*_1.fastq.gz' ! -iname '*_2.fastq.gz' \
      ! -iname '*_R1*.fastq.gz' ! -iname '*_R2*.fastq.gz' \
    | while read -r f; do
      echo -e "${study}\tAMBIGUOUS_BARE\t$(basename "$f")" >> "$REPORT"
      log "  AMBIGUOUS: $(basename "$f") has no _1/_2/R1/R2 marker in a PAIRED study - needs manual check"
    done
  fi
done

echo
log "=== NORMALIZATION PASS COMPLETE ==="
echo
echo "--- Issues requiring manual attention ---"
if [[ $(wc -l < "$REPORT") -le 1 ]]; then
  echo "  none - clean pass"
else
  column -t -s$'\t' "$REPORT" 2>/dev/null || cat "$REPORT"
fi
echo
echo "Full report: $REPORT"
echo "Next: resolve anything listed above, then run 01_build_manifests.sh"
