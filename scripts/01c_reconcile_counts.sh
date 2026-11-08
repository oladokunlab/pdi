#!/bin/bash
# =============================================================================
# 01c_reconcile_counts.sh - find WHY a study's file count disagrees with metadata
#
# HOW TO RUN: login node, read-only. Changes nothing, deletes nothing.
#   bash 01c_reconcile_counts.sh                 # all studies
#   bash 01c_reconcile_counts.sh Ma_2024         # one study, verbose
#
# WHY: 01_build_manifests only tells me a study is INCOMPLETE, not why. Getting
# this wrong in either direction is expensive - extra files mean samples in my
# feature table that are not in the metadata (QIIME 2 will error, or worse,
# get silently dropped at merge), and missing files mean a study is
# under-represented. This tells me which it is and where the extras come from.
#
# Expected files per study:
#   paired         -> n_runs x 2   (_1 and _2)
#   single_merged  -> n_runs x 1
#   ccs            -> n_runs x 1
# =============================================================================

source "$(dirname "$0")/00_config.sh"

if [[ -z "${LOG_DIR:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh did not load correctly." >&2
  exit 1
fi

TARGET="${1:-}"

printf '%-24s %8s %8s %7s  %-14s %s\n' study actual expect diff layout note
printf '%.0s-' {1..95}; echo

tail -n +2 "$STUDY_TABLE" | while IFS=$'\t' read -r study n_runs region platform layout pf pr status; do
  [[ -z "$study" ]] && continue
  [[ -n "$TARGET" && "$study" != "$TARGET" ]] && continue

  sdir="${RAW_DIR}/${study}"
  if [[ ! -d "$sdir" ]]; then
    printf '%-24s %8s %8s %7s  %-14s %s\n' "$study" 0 "?" "" "$layout" "NO DIRECTORY"
    continue
  fi

  actual=$(find "$sdir" -iname '*.fastq.gz' -o -iname '*.fastq' 2>/dev/null | wc -l)
  case "$layout" in
    paired)  expect=$((n_runs * 2)) ;;
    *)       expect=$n_runs ;;
  esac
  diff=$((actual - expect))

  note=""
  if [[ "$diff" -ne 0 && "$n_runs" -gt 0 ]]; then
    # integer files-per-run ratio is a strong hint at the cause
    ratio=$(awk -v a="$actual" -v n="$n_runs" 'BEGIN{printf "%.2f", a/n}')
    note="${ratio} files/run"
  fi

  printf '%-24s %8s %8s %+7s  %-14s %s\n' "$study" "$actual" "$expect" "$diff" "$layout" "$note"

  # ---- verbose breakdown, only when a single study was requested -----------
  if [[ -n "$TARGET" ]]; then
    echo
    echo "  --- distinct run accessions (mate markers stripped) ---"
    find "$sdir" -iname '*.fastq*' -printf '%f\n' \
      | sed -E 's/_[Rr]?[12]\.fastq.*$//; s/\.fastq.*$//; s/\.sra$//' \
      | sort -u > /tmp/_runs_$$.txt
    echo "  distinct accessions on disk: $(wc -l < /tmp/_runs_$$.txt)  (metadata says ${n_runs})"
    echo
    echo "  --- accessions on disk but NOT in metadata (candidates to remove) ---"
    awk -F'\t' -v s="$study" 'NR>2 && $2==s {print $1}' "$METADATA" | sort -u > /tmp/_meta_$$.txt
    comm -23 /tmp/_runs_$$.txt /tmp/_meta_$$.txt | head -40
    n_extra=$(comm -23 /tmp/_runs_$$.txt /tmp/_meta_$$.txt | wc -l)
    echo "  ($n_extra extra)"
    echo
    echo "  --- in metadata but MISSING from disk ---"
    comm -13 /tmp/_runs_$$.txt /tmp/_meta_$$.txt | head -40
    n_miss=$(comm -13 /tmp/_runs_$$.txt /tmp/_meta_$$.txt | wc -l)
    echo "  ($n_miss missing)"
    echo
    echo "  --- files per accession (top offenders) ---"
    find "$sdir" -iname '*.fastq*' -printf '%f\n' \
      | sed -E 's/_[Rr]?[12]\.fastq.*$//; s/\.fastq.*$//; s/\.sra$//' \
      | sort | uniq -c | sort -rn | head -10
    rm -f /tmp/_runs_$$.txt /tmp/_meta_$$.txt
    echo
  fi
done

echo
echo "Studies with diff != 0 need resolving BEFORE stage 02."
echo "Re-run against one study for a full breakdown, e.g.:"
echo "  bash $(basename "$0") Ma_2024"
