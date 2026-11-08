#!/bin/bash
# =============================================================================
# 01e_reconcile_accessions.sh - reconcile files on disk against the metadata
#
# HOW TO RUN: login node. REPORT MODE BY DEFAULT - nothing moves.
#   bash 01e_reconcile_accessions.sh                 # report only
#   bash 01e_reconcile_accessions.sh --quarantine    # move extras out
#
# WHY: file COUNTS matching is not the same as the right FILES being present.
# Melaku_2024 had the correct number of runs on disk and not a single one of
# them was a run the metadata actually calls for - the download had pulled a
# different slice of the BioProject entirely. A count-based check said
# "8 files, fine". An accession-based check says "zero of these are yours".
# That is the difference this script exists to catch.
#
# Two failure modes it separates:
#   EXTRA   - on disk, not in the metadata. Typically occurs when the whole
#             BioProject was downloaded but only some arms are in scope. These
#             are quarantined, never deleted.
#   MISSING - in the metadata, not on disk. These must be downloaded before the
#             study can be used. The script writes a ready-to-run accession
#             list per study.
#
# NOTHING IS EVER DELETED. Extras move to raw_data_excluded/<study>/ so the
# decision is reversible and auditable.
# =============================================================================

source "$(dirname "$0")/00_config.sh"

if [[ -z "${LOG_DIR:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh did not load correctly." >&2
  exit 1
fi

QUARANTINE=0
[[ "${1:-}" == "--quarantine" ]] && QUARANTINE=1

EXCL_DIR="${PDI_ROOT}/raw_data_excluded"
REDL_DIR="${PDI_ROOT}/redownload_lists"
mkdir -p "$REDL_DIR"

if [[ "$QUARANTINE" -eq 1 ]]; then
  log "QUARANTINE MODE - extras will be MOVED to ${EXCL_DIR}"
else
  log "REPORT MODE - nothing will be moved. Add --quarantine to act."
fi
echo

printf '%-24s %6s %6s %6s %6s  %s\n' STUDY DISK META EXTRA MISS STATUS
printf '%.0s-' {1..80}; echo

tot_extra=0; tot_missing=0; studies_clean=0; studies_bad=0

while IFS=$'\t' read -r study n_runs region platform layout pf pr status pp orient; do
  [[ -z "$study" ]] && continue
  sdir="${RAW_DIR}/${study}"
  [[ -d "$sdir" ]] || { printf '%-24s %6s %6s %6s %6s  %s\n' "$study" 0 "?" "-" "-" "NO_DIRECTORY"; continue; }

  # accessions present on disk. Strip mate markers AND the .sra infix that
  # fastq-dump leaves in filenames, or every file looks like an extra.
  find "$sdir" -maxdepth 1 -iname '*.fastq*' -printf '%f\n' 2>/dev/null \
    | sed -E 's/_[Rr]?[12]\.fastq.*$//; s/\.fastq.*$//; s/\.sra$//' \
    | sort -u > /tmp/_disk_$$.txt

  # accessions the metadata says belong to this study
  awk -F'\t' -v s="$study" 'NR>2 && $2==s {print $1}' "$METADATA" \
    | tr -d '\r' | sort -u > /tmp/_meta_$$.txt

  n_disk=$(wc -l < /tmp/_disk_$$.txt)
  n_meta=$(wc -l < /tmp/_meta_$$.txt)

  comm -23 /tmp/_disk_$$.txt /tmp/_meta_$$.txt > /tmp/_extra_$$.txt
  comm -13 /tmp/_disk_$$.txt /tmp/_meta_$$.txt > /tmp/_miss_$$.txt
  n_extra=$(wc -l < /tmp/_extra_$$.txt)
  n_miss=$(wc -l < /tmp/_miss_$$.txt)

  if   [[ "$n_extra" -eq 0 && "$n_miss" -eq 0 ]]; then verdict="CLEAN"; studies_clean=$((studies_clean+1))
  elif [[ "$n_miss" -gt 0 && "$n_extra" -gt 0 ]]; then verdict="WRONG_SUBSET_DOWNLOADED"; studies_bad=$((studies_bad+1))
  elif [[ "$n_miss" -gt 0 ]];                     then verdict="INCOMPLETE_needs_download"; studies_bad=$((studies_bad+1))
  else                                                 verdict="EXTRAS_only"; studies_bad=$((studies_bad+1))
  fi

  printf '%-24s %6s %6s %6s %6s  %s\n' "$study" "$n_disk" "$n_meta" "$n_extra" "$n_miss" "$verdict"

  tot_extra=$((tot_extra + n_extra))
  tot_missing=$((tot_missing + n_miss))

  # write a re-download list for anything missing
  if [[ "$n_miss" -gt 0 ]]; then
    cp /tmp/_miss_$$.txt "${REDL_DIR}/${study}_missing.txt"
  fi

  # quarantine extras
  if [[ "$n_extra" -gt 0 && "$QUARANTINE" -eq 1 ]]; then
    mkdir -p "${EXCL_DIR}/${study}"
    while read -r acc; do
      [[ -z "$acc" ]] && continue
      for f in "$sdir/$acc"*; do
        [[ -e "$f" ]] && mv -n "$f" "${EXCL_DIR}/${study}/"
      done
    done < /tmp/_extra_$$.txt
    log "  ${study}: quarantined ${n_extra} accessions -> ${EXCL_DIR}/${study}/"
  fi

  rm -f /tmp/_disk_$$.txt /tmp/_meta_$$.txt /tmp/_extra_$$.txt /tmp/_miss_$$.txt

done < <(tail -n +2 "$STUDY_TABLE" | tr -d '\r')

echo
log "clean studies: ${studies_clean}   studies needing attention: ${studies_bad}"
log "total extra accessions: ${tot_extra}   total missing: ${tot_missing}"
echo
if [[ "$tot_missing" -gt 0 ]]; then
  echo "Re-download lists written to ${REDL_DIR}/ :"
  ls -1 "$REDL_DIR" 2>/dev/null | sed 's/^/  /'
  echo
  echo "To fetch the missing runs for a study (LOGIN NODE):"
  echo "  cd ${RAW_DIR}/<STUDY>"
  echo "  nohup bash -c 'while read a; do fasterq-dump --split-files --threads 8 \"\$a\" -O .; gzip \"\$a\"_1.fastq \"\$a\"_2.fastq; done < ${REDL_DIR}/<STUDY>_missing.txt' > ${LOG_DIR}/<STUDY>_redl.log 2>&1 &"
fi
[[ "$QUARANTINE" -eq 0 && "$tot_extra" -gt 0 ]] && \
  echo && echo "Re-run with --quarantine to move the ${tot_extra} extra accessions out of raw_data/."
exit 0
