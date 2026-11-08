#!/usr/bin/env bash
# =============================================================================
# verify_dependencies.sh
#
# Checks every key output file against everything it depends on. An output
# older than one of its own inputs means it was built before that input's
# last change - i.e. it is stale and needs rerunning.
#
# RUN:  bash verify_dependencies.sh
#
# This does not re-run anything or modify any file. It only compares
# modification timestamps and reports PASS/FAIL/MISSING for every link in the
# dependency chain, so nothing has to be taken on trust or reconstructed from
# memory across a long session.
# =============================================================================
set -uo pipefail

ROOT="/scratch/user/$USER/pdi"
SCR="${ROOT}/scripts"
GHI="${ROOT}/r_analysis/ghi"
GDI="${ROOT}/r_analysis/gdi"
GEN="${ROOT}/r_analysis/genus"

PASS=0; FAIL=0; MISS=0

mtime() { stat -c %Y "$1" 2>/dev/null; }
fmt()   { stat -c '%y' "$1" 2>/dev/null | cut -d. -f1; }

check() {
  local label="$1"; shift
  local output="$1"; shift
  local inputs=("$@")

  if [ ! -e "$output" ]; then
    echo "  MISSING  $label"
    echo "           output not found: $output"
    MISS=$((MISS+1)); return
  fi
  local out_t; out_t=$(mtime "$output")
  local worst_in=""; local worst_t=0; local any_missing=0

  for inp in "${inputs[@]}"; do
    if [ ! -e "$inp" ]; then
      echo "  MISSING  $label"
      echo "           input not found: $inp"
      any_missing=1
      continue
    fi
    local in_t; in_t=$(mtime "$inp")
    if [ "$in_t" -gt "$worst_t" ]; then worst_t=$in_t; worst_in=$inp; fi
  done
  if [ "$any_missing" -eq 1 ]; then MISS=$((MISS+1)); return; fi

  if [ "$out_t" -ge "$worst_t" ]; then
    echo "  PASS     $label"
    PASS=$((PASS+1))
  else
    echo "  FAIL     $label"
    echo "           output : $(fmt "$output")   $output"
    echo "           stale vs: $(fmt "$worst_in")   $worst_in"
    FAIL=$((FAIL+1))
  fi
}

echo "============================================================================"
echo "DEPENDENCY VERIFICATION"
echo "============================================================================"
echo ""
echo "Level 0 - corrected metadata (both classification corrections applied)"
check "Metadata is the final corrected version" \
  "${SCR}/GDI_Metadata.xlsx" \
  "${SCR}/GDI_Metadata_backup_20260819_152409.xlsx"

echo ""
echo "Level 1 - genus agglomeration"
check "genus_counts.tsv built after final metadata" \
  "${GEN}/genus_counts.tsv" \
  "${SCR}/GDI_Metadata.xlsx"

echo ""
echo "Level 2 - batch correction (age covariate)"
check "genus_clr_batch_adjusted.tsv built after genus_counts.tsv and metadata" \
  "${GEN}/genus_clr_batch_adjusted.tsv" \
  "${GEN}/genus_counts.tsv" "${SCR}/GDI_Metadata.xlsx"

echo ""
echo "Level 3 - GHI construction (morphology fix, Performance split, new domains)"
check "ghi_scores_pca.tsv built after metadata (biomarker sheet)" \
  "${GHI}/ghi_scores_pca.tsv" \
  "${SCR}/GDI_Metadata.xlsx"

echo ""
echo "Level 4 - GHI validation (must postdate the 04e domain-list sync fix)"
check "mixed_model_validation.txt built after ghi_scores_pca.tsv" \
  "${GHI}/mixed_model_validation.txt" \
  "${GHI}/ghi_scores_pca.tsv"

echo ""
echo "Level 5 - GDI core chain (all must postdate CLR matrix and GHI scores)"
check "maaslin2_all_results.tsv built after CLR matrix and GHI scores" \
  "${GDI}/maaslin2_all_results.tsv" \
  "${GEN}/genus_clr_batch_adjusted.tsv" "${GHI}/ghi_scores_pca.tsv"

check "gdi_model.rds (LASSO) built after MaAsLin2 output" \
  "${GDI}/gdi_model.rds" \
  "${GDI}/maaslin2_all_results.tsv"

check "gdi_stability.tsv built no earlier than gdi_model.rds" \
  "${GDI}/gdi_stability.tsv" \
  "${GDI}/gdi_model.rds"

check "gdi_scores.tsv (GDI construction) built after gdi_stability.tsv" \
  "${GDI}/gdi_scores.tsv" \
  "${GDI}/gdi_stability.tsv"

check "loso_summary.txt built after gdi_scores.tsv" \
  "${GDI}/loso_summary.txt" \
  "${GDI}/gdi_scores.tsv"

echo ""
echo "Level 6 - sensitivity analyses (must all postdate gdi_model.rds, i.e. the"
echo "          corrected classification + age covariate + new domains)"
for f in \
  "sensitivity_analyses/stratified_test.txt|Stratified arm-subset test (06e)" \
  "sensitivity_analyses/gdi_by_class.txt|GDI fitted per class (14)" \
  "sensitivity_analyses/three_way_axis.txt|Three-way unsupervised axis (16)" \
  "control_vs_treated.txt|Control vs everything pooled (20)" \
  "target_comparison_domains.txt|Four-target comparison (21)" \
  "taxon_metaanalysis.txt|Per-taxon meta-analysis (18)"
do
  path="${f%%|*}"; label="${f##*|}"
  check "$label" "${GDI}/${path}" "${GDI}/gdi_model.rds"
done

echo ""
echo "Level 7 - sparse PLS (fallback, spls package)"
check "spls_fallback summary built after GHI scores and CLR matrix" \
  "${GDI}/spls_fallback/spls_summary.txt" \
  "${GHI}/ghi_scores_pca.tsv" "${GEN}/genus_clr_batch_adjusted.tsv"

echo ""
echo "============================================================================"
echo "SUMMARY:  PASS=${PASS}   FAIL=${FAIL}   MISSING=${MISS}"
echo "============================================================================"
if [ "$FAIL" -gt 0 ]; then
  echo ""
  echo "Any FAIL above means that file was built BEFORE one of its own inputs last"
  echo "changed - it does not reflect the current, corrected data and must be"
  echo "regenerated before being cited in the report."
fi
if [ "$MISS" -gt 0 ]; then
  echo ""
  echo "Any MISSING above means either the output was never produced, or one of"
  echo "its expected inputs is not where this script looked - check the exact path."
fi
