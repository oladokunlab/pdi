#!/bin/bash
# =============================================================================
# 00_preflight_check.sh
#
# Validates the pipeline configuration and input files before submitting any
# compute jobs. Run on the login node; completes in seconds.
#
# Usage:
#   bash 00_preflight_check.sh
# =============================================================================

source "$(dirname "$0")/00_config.sh"

FAIL=0
ok()   { echo "  [ OK ] $*"; }
bad()  { echo "  [FAIL] $*"; FAIL=1; }
warn() { echo "  [WARN] $*"; }

echo "=========================================="
echo " PDI PIPELINE PREFLIGHT"
echo "=========================================="

echo
echo "-- Environment --"
[[ -x "$QIIME"      ]] && ok "qiime at $QIIME"        || bad "qiime NOT FOUND at $QIIME"
[[ -s "$CLASSIFIER" ]] && ok "SILVA classifier found" || bad "classifier NOT FOUND at $CLASSIFIER"
command -v biom >/dev/null && ok "biom on PATH"       || warn "biom not on PATH (needed in stage 04)"

echo
echo "-- Metadata --"
if [[ -s "$METADATA" ]]; then
  n=$(( $(wc -l < "$METADATA") ))
  # Expected line count is derived from studies.tsv rather than hardcoded, so
  # this check stays valid as the corpus size changes.
  exp_samples=$(tail -n +2 "$STUDY_TABLE" | tr -d '\r' | awk -F'\t' '{s+=$2} END{print s}')
  exp_lines=$((exp_samples + 2))
  if [[ "$n" -eq "$exp_lines" ]]; then
    ok "metadata: ${n} lines = header + #q2:types + ${exp_samples} samples"
  else
    bad "metadata has ${n} lines, expected ${exp_lines} (${exp_samples} samples per studies.tsv)"
  fi
  head -2 "$METADATA" | tail -1 | grep -q '#q2:types' \
    && ok "#q2:types row present on line 2" \
    || bad "#q2:types row MISSING from line 2"
  # QIIME 2 requires LF line endings; CRLF causes the final column name to be
  # parsed with a trailing carriage return.
  if grep -qU $'\r' "$METADATA"; then
    warn "CR characters present. Fix with: sed -i 's/\\r\$//' \"\$METADATA\""
  else
    ok "line endings clean (LF only)"
  fi
else
  bad "metadata NOT FOUND at $METADATA"
fi

echo
echo "-- Study table --"
if [[ -s "$STUDY_TABLE" ]]; then
  n=$(( $(wc -l < "$STUDY_TABLE") - 1 ))
  ok "studies.tsv found (${n} studies)"
  total=$(tail -n +2 "$STUDY_TABLE" | tr -d '\r' | awk -F'\t' '{s+=$2} END{print s}')
  ok "run counts sum to ${total} across ${n} studies"

  # Every study in studies.tsv must also appear in the metadata file, and
  # vice versa - a mismatch here causes the merge in stage 03 to fail.
  if [[ -s "$METADATA" ]]; then
    tail -n +2 "$STUDY_TABLE" | tr -d '\r' | cut -f1 | sort -u > /tmp/_st_$$.txt
    awk -F'\t' 'NR>2{gsub(/\r/,""); print $2}' "$METADATA" | sort -u > /tmp/_mt_$$.txt
    only_st=$(comm -23 /tmp/_st_$$.txt /tmp/_mt_$$.txt | tr '\n' ' ')
    only_mt=$(comm -13 /tmp/_st_$$.txt /tmp/_mt_$$.txt | tr '\n' ' ')
    [[ -z "$only_st" ]] && ok "every studies.tsv study present in metadata" \
                        || bad "in studies.tsv but NOT metadata: ${only_st}"
    [[ -z "$only_mt" ]] && ok "every metadata study present in studies.tsv" \
                        || bad "in metadata but NOT studies.tsv: ${only_mt}"
    rm -f /tmp/_st_$$.txt /tmp/_mt_$$.txt
  fi
  # The SLURM array size in stage 02 must match the study count; otherwise
  # array indices beyond the study count exit silently without processing
  # anything.
  arr=$(grep -oP '(?<=--array=1-)\d+' "$(dirname "$0")/02_denoise_per_study.sh" || echo 0)
  [[ "$arr" -eq "$n" ]] && ok "SBATCH array size (${arr}) matches study count" \
                        || bad "SBATCH array is 1-${arr} but there are ${n} studies"
else
  bad "studies.tsv NOT FOUND at $STUDY_TABLE"
fi

echo
echo "-- Raw data --"
missing=0; present=0
while IFS=$'\t' read -r study n region platform layout pf pr status; do
  d="${RAW_DIR}/${study}"
  if [[ -d "$d" ]]; then
    c=$(find "$d" -maxdepth 1 \( -iname '*.fastq' -o -iname '*.fastq.gz' \) 2>/dev/null | wc -l)
    if [[ "$c" -eq 0 ]]; then
      echo "    ${study}: directory exists but 0 FASTQ"; missing=$((missing+1))
    else
      present=$((present+1))
    fi
  else
    echo "    ${study}: NO DIRECTORY (${status})"; missing=$((missing+1))
  fi
done < <(tail -n +2 "$STUDY_TABLE")
ok "${present} studies have FASTQ on disk"
[[ "$missing" -gt 0 ]] && warn "${missing} studies missing data - check status column for reason"

echo
echo "-- Studies on hold --"
tail -n +2 "$STUDY_TABLE" | awk -F'\t' '$8 ~ /^HOLD/ {printf "    %-22s %-4s %s\n", $1, $2, $8}'
tail -n +2 "$STUDY_TABLE" | awk -F'\t' '$8 ~ /^WARN/ {printf "    %-22s %-4s %s\n", $1, $2, $8}'

echo
echo "=========================================="
[[ "$FAIL" -eq 0 ]] && echo " PREFLIGHT PASSED - safe to run stage 01" \
                    || echo " PREFLIGHT FAILED - fix the [FAIL] items above"
echo "=========================================="
exit $FAIL
