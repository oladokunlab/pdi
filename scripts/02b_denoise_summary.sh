#!/bin/bash
# =============================================================================
# 02b_denoise_summary.sh - pull DADA2 retention for every study into one table
#
# HOW TO RUN: login node (cheap - just exports .qza files).
#   bash 02b_denoise_summary.sh | tee ../logs/02b_denoise_summary.log
#
# WHY: opening 51 separate denoising-stats.qzv files in a browser is not a
# realistic quality check. This flattens all of them into one table to
# see at a glance which studies lost reads and where.
#
# The columns DADA2 reports, in order, are:
#   input -> filtered -> denoised -> merged -> non-chimeric
# and the place a study bleeds reads tells me WHY:
#   big loss at 'filtered'      -> quality worse than max-EE 2.0 tolerates
#   big loss at 'merged'        -> R1/R2 do not overlap (amplicon too long for
#                                  the read length) - this is what killed
#                                  Wassie_2021
#   big loss at 'non-chimeric'  -> over-amplification during PCR
#
# Typical benchmark range: approximately 87-96% retained at filtering,
# 86-95% non-chimeric. Studies under 50% overall retention are flagged for
# closer review.
#
# Particular attention is warranted for two groups:
#   - the 8 studies where cutadapt was skipped (primers pre-removed by the
#     submitter) - a misclassification here causes retention to collapse
#   - the 30 studies on NovaSeq/NextSeq-class platforms with 4-level binned
#     quality scores, where DADA2's error model has far less information to
#     work with than it was designed for
# =============================================================================

source "$(dirname "$0")/00_config.sh"
if [[ -z "${LOG_DIR:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh did not load correctly." >&2; exit 1
fi
activate_qiime

OUT="${LOG_DIR}/denoise_stats_all.tsv"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

printf 'study\tn_samples\tinput\tfiltered\tdenoised\tmerged\tnon_chimeric\tpct_filtered\tpct_merged\tpct_final\n' > "$OUT"

while IFS=$'\t' read -r study n_runs region platform layout pf pr status pp orient; do
  [[ -z "$study" ]] && continue
  stats="${WORK_DIR}/${study}/denoising-stats.qza"
  [[ -s "$stats" ]] || { echo "  MISSING stats: ${study}" >&2; continue; }

  rm -rf "${TMP}/x"; mkdir -p "${TMP}/x"
  $QIIME tools export --input-path "$stats" --output-path "${TMP}/x" >/dev/null 2>&1 || {
    echo "  EXPORT FAILED: ${study}" >&2; continue; }

  f="${TMP}/x/stats.tsv"
  [[ -s "$f" ]] || continue

  # stats.tsv has a header row then a '#q2:types' row, then one row per sample.
  # Column names differ between denoise-paired (has 'merged') and
  # denoise-single / denoise-ccs (no merge step), so columns are looked up by name
  # rather than by position.
  awk -F'\t' -v S="$study" '
    NR==1 {
      for (i=1;i<=NF;i++) h[$i]=i
      ci=h["input"]; cf=h["filtered"]; cd=h["denoised"]
      cm=("merged" in h) ? h["merged"] : 0
      cn=("non-chimeric" in h) ? h["non-chimeric"] : 0
      next
    }
    $1 ~ /^#/ { next }
    {
      n++
      i+=$ci; fi+=$cf; d+=$cd
      if (cm) m+=$cm; if (cn) nc+=$cn
    }
    END {
      if (n==0) exit
      mm = (cm ? m : d)          # no merge step for single/ccs
      nn = (cn ? nc : mm)
      printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\n",
             S, n, i, fi, d, mm, nn,
             (i? 100*fi/i:0), (i? 100*mm/i:0), (i? 100*nn/i:0)
    }' "$f" >> "$OUT"
done < <(tail -n +2 "$STUDY_TABLE" | tr -d '\r')

echo
column -t "$OUT" 2>/dev/null || cat "$OUT"

echo
echo "=============================================="
echo " STUDIES NEEDING A LOOK (final retention <50%)"
echo "=============================================="
awk -F'\t' 'NR>1 && $10+0 < 50 {printf "  %-24s final=%.1f%%  (filtered=%.1f%% merged=%.1f%%)\n", $1,$10,$8,$9}' "$OUT"
n_bad=$(awk -F'\t' 'NR>1 && $10+0 < 50' "$OUT" | wc -l)
[[ "$n_bad" -eq 0 ]] && echo "  none - all studies above 50%"

echo
echo "corpus totals:"
awk -F'\t' 'NR>1{s+=$2; i+=$3; f+=$7} END{printf "  %d samples, %d reads in, %d retained (%.1f%%)\n", s, i, f, (i?100*f/i:0)}' "$OUT"
echo
echo "full table: $OUT"
