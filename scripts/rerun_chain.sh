#!/bin/bash
# Re-run the analysis chain in dependency order so every output postdates its inputs.
set -u
export PDI_ROOT="${PDI_ROOT:-/scratch/user/$USER/pdi}"
export R_LIBS_USER="${R_LIBS_USER:-/scratch/user/$USER/R/library}"
RS=/sw/eb/sw/R/4.5.2-gfbf-2025b/bin/Rscript
cd "$PDI_ROOT/scripts/R_scripts" || exit 1
STAMP=$(date +%Y%m%d_%H%M)
LOG="$PDI_ROOT/qc_reports/rerun_chain_${STAMP}.log"
mkdir -p "$PDI_ROOT/qc_reports"
MARK="$PDI_ROOT/qc_reports/.chain_start_${STAMP}"; touch "$MARK"

# Numbered scripts only, version-sorted. 03 excluded: the genus table is final
# and MMUPHin need not re-run. 19 excluded: mixOmics is not installable here.
LIST=$(ls [0-9][0-9]*.R 2>/dev/null | grep -v '^03' | grep -v '^19_sparse_pls.R' | LC_ALL=C sort)

PASS=0; FAIL=0; FAILED=""
echo "chain start $(date)" | tee "$LOG"
for f in $LIST; do
  t0=$SECONDS
  if $RS "$f" >>"$LOG" 2>&1; then
    PASS=$((PASS+1)); printf '  ok   %-42s %4ds\n' "$f" $((SECONDS-t0)) | tee -a "$LOG"
  else
    FAIL=$((FAIL+1)); FAILED="$FAILED $f"; printf '  FAIL %-42s %4ds\n' "$f" $((SECONDS-t0)) | tee -a "$LOG"
  fi
done

echo "" | tee -a "$LOG"
echo "passed $PASS, failed $FAIL" | tee -a "$LOG"
[ -n "$FAILED" ] && echo "failed:$FAILED" | tee -a "$LOG"

echo "" | tee -a "$LOG"
# 05_maaslin2.R names its output directory after the covariate set it is given.
# The loop above runs the default; the other three configurations are run here
# so that every deposited MaAsLin2 directory comes from this pass.
for cfg in "no no no" "yes no no" "no yes yes"; do
  set -- $cfg
  if AGE_COVARIATE=$1 CLASS_COVARIATE=$2 SEGMENT_COVARIATE=$3 $RS 05_maaslin2.R >>"$LOG" 2>&1; then
    printf '  ok   05_maaslin2.R age=%s class=%s seg=%s\n' "$1" "$2" "$3" | tee -a "$LOG"
  else
    printf '  FAIL 05_maaslin2.R age=%s class=%s seg=%s\n' "$1" "$2" "$3" | tee -a "$LOG"
  fi
done

echo "--- outputs NOT refreshed by this pass ---" | tee -a "$LOG"
find "$PDI_ROOT/r_analysis" -type f \( -name '*.txt' -o -name '*.tsv' \) \
  ! -newer "$MARK" -printf '%TY-%Tm-%Td_%TH:%TM  %p\n' | sort | tee -a "$LOG"
find "$PDI_ROOT/r_analysis" -type f \( -name '*.txt' -o -name '*.tsv' \) \
  ! -newer "$MARK" | wc -l | xargs echo "stale count:" | tee -a "$LOG"
mv -f "$MARK" "$PDI_ROOT/qc_reports/.last_chain_start"
echo "log: $LOG"
