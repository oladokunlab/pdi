#!/usr/bin/env bash
# =============================================================================
# cleanup_nonessential.sh
#
# Moves non-essential scripts and stale outputs into _removed/ so the working
# directory reflects only what the project actually depends on. Nothing is
# deleted - everything stays recoverable in _removed/ with a manifest saying
# why each item was set aside.
#
# RUN:  bash cleanup_nonessential.sh          (dry run - lists only)
#       bash cleanup_nonessential.sh --apply  (performs the moves)
#
# -----------------------------------------------------------------------------
# HOW THINGS WERE CLASSIFIED
# -----------------------------------------------------------------------------
# ESSENTIAL   - a reported number depends on it, or it is required to
#               reproduce one. Never moved.
# MID         - not needed to reproduce a number, but documents a decision,
#               a correction, or a negative result that appears in the write-up.
#               Never moved.
# NON         - superseded by a later script, a dead end nothing cites, or a
#               one-off diagnostic whose finding is already recorded elsewhere.
#               Moved to _removed/.
#
# A script that produced a NULL result is NOT non-essential: the nulls from
# binary encoding, percentile normalisation, the elastic-net sweep and random
# forest are all reported findings, so those scripts stay.
# =============================================================================
set -uo pipefail

ROOT="/scratch/user/$USER/pdi"
R="${ROOT}/scripts/R_scripts"
S="${ROOT}/scripts"
GDI="${ROOT}/r_analysis/gdi"
STAMP=$(date +%Y%m%d_%H%M)
REM="${ROOT}/_removed_${STAMP}"

APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1

mkdir -p "${REM}/R_scripts" "${REM}/bash" "${REM}/outputs" 2>/dev/null

MANIFEST="${REM}/WHY_THESE_WERE_REMOVED.txt"
n_moved=0; n_missing=0

move_it() {
  local src="$1" dest="$2" why="$3"
  if [[ ! -e "$src" ]]; then
    echo "  (not present) $(basename "$src")"
    n_missing=$((n_missing+1)); return
  fi
  echo "  MOVE  $(basename "$src")"
  echo "        reason: $why"
  n_moved=$((n_moved+1))
  if [[ $APPLY -eq 1 ]]; then
    mv "$src" "$dest"/ 2>/dev/null
    printf '%-40s %s\n' "$(basename "$src")" "$why" >> "$MANIFEST"
  fi
}

echo "============================================================"
echo "SUPERSEDED BY A LATER SCRIPT"
echo "============================================================"
move_it "${R}/06b_target_comparison.R"   "${REM}/R_scripts" "superseded by 21_target_comparison_domains.R"
move_it "${R}/06c_lambda_comparison.R"   "${REM}/R_scripts" "superseded by 24_elastic_net.R (full alpha sweep)"
move_it "${R}/06d_panel_size_test.R"     "${REM}/R_scripts" "superseded by 21_target_comparison_domains.R"
move_it "${R}/10_unsupervised_gdi.R"     "${REM}/R_scripts" "superseded by 12_gdi_three_ways.R"
move_it "${R}/11_axis_loso.R"            "${REM}/R_scripts" "superseded by 16_three_way_axis.R"
move_it "${R}/fix_labels.R"              "${REM}/R_scripts" "superseded by fix_classes.R and fix_classes2.R"

echo ""
echo "============================================================"
echo "DEAD ENDS - no reported result depends on these"
echo "============================================================"
move_it "${S}/05_diversity.sh"                  "${REM}/bash"      "never produced output; alpha diversity not used in any result"
move_it "${R}/09_absolute_discriminant_test.R"  "${REM}/R_scripts" "single-sample deployability test; not cited in the report"

echo ""
echo "============================================================"
echo "ONE-OFF DIAGNOSTICS - finding already recorded elsewhere"
echo "============================================================"
move_it "${R}/check_dropped_studies.R"    "${REM}/R_scripts" "one-off: traced 40->40 study count, no drop found"
move_it "${R}/check_44_to_40.R"           "${REM}/R_scripts" "one-off: identified 3 covariate-dropped studies, now fixed"
move_it "${R}/check_44_to_43.R"           "${REM}/R_scripts" "one-off: identified Liu_2026 lost at zero-replacement"
move_it "${R}/check_why_dropped.R"        "${REM}/R_scripts" "one-off: identified which field each dropped study was missing"
move_it "${R}/check_ma2024_rows.R"        "${REM}/R_scripts" "one-off: exposed the Ma_2024 two-experiment conflation"
move_it "${R}/check_ma2024_sequencing.R"  "${REM}/R_scripts" "one-off: confirmed which Ma_2024 trial had sequencing"
move_it "${R}/check_header_raw.R"         "${REM}/R_scripts" "one-off: diagnosed the duplicated #q2:types header line"
move_it "${R}/check_tsv_integrity.R"      "${REM}/R_scripts" "one-off: post-corruption integrity check, issue resolved"
move_it "${R}/check_age_plausibility.R"   "${REM}/R_scripts" "superseded by verify_ages_vs_sources.R"
move_it "${R}/audit_metadata.R"           "${REM}/R_scripts" "superseded by verify_metadata_full.R"
move_it "${R}/audit_remaining_fixes.R"    "${REM}/R_scripts" "one-off: scoped the age fixes, all now applied"
move_it "${R}/show_arm_mismatches.R"      "${REM}/R_scripts" "one-off: scoped the arm-name fixes, all now applied"

echo ""
echo "============================================================"
echo "STALE OUTPUT"
echo "============================================================"
move_it "${GDI}/gdi_scores.tsv" "${REM}/outputs" "Aug 18 file; superseded - no panel is constructed on the final data"

echo ""
echo "============================================================"
echo "KEPT (for reference - these are NOT moved)"
echo "============================================================"
cat <<'KEPT'
  ESSENTIAL - a reported number depends on these:
    00_config.sh, 00_r_config.R, 01_build_manifests.sh, 02_denoise_per_study.sh
    03_merge_and_taxonomy.sh, 03a_merge_and_pcoa.sh, 03b_build_phylogeny.sh
    03c_refresh_pcoa_permanova.sh, 04_filter_and_export.sh
    01_genus_agglomeration.R, 03_batch_correction.R, 04b_ghi_pca.R
    04e_mixed_model_validation.R, 04g_domain_correlation.R
    05_maaslin2.R, 06_lasso_selection.R, 07_gdi_construct.R
    08_loso_validate.R, make_figures.R

  MID - documents a correction, decision, or reported negative result:
    fix_classes.R, fix_classes2.R, fix_missing_ages.R, fix_ages_final.R
    fix_ma2024_split.R, rename_ma2024_back.R, fix_arm_names.R
    12_gdi_three_ways.R, 13_ghi_permanova.R, 14_gdi_by_class.R
    16_three_way_axis.R, 17_maaslin_class.R, 18_taxon_metaanalysis.R
    19_sparse_pls.R (attempted), 19b_sparse_pls_fallback.R
    20_control_vs_treated.R, 21_target_comparison_domains.R
    22_binary_presence_lasso.R, 23_percentile_normalization.R
    24_elastic_net.R, 25_random_forest.R
    06e_stratified_test.R, 06f_challenge_gdi.R
    variance_partition.R, units_audit.R, barrier_and_coverage.R
    verify_metadata_full.R, verify_dependencies.sh
    verify_ages_vs_sources.R
    every GDI_Metadata*_backup_*.xlsx and *_backup_*.tsv  (audit trail)
KEPT

echo ""
echo "============================================================"
if [[ $APPLY -eq 1 ]]; then
  echo "MOVED: ${n_moved} item(s) -> ${REM}"
  echo "Manifest: ${MANIFEST}"
else
  rmdir "${REM}/R_scripts" "${REM}/bash" "${REM}/outputs" "${REM}" 2>/dev/null
  echo "DRY RUN: ${n_moved} item(s) would move, ${n_missing} not present."
  echo "Re-run with --apply to perform the moves."
fi
echo "============================================================"
