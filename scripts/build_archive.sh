#!/usr/bin/env bash
# =============================================================================
# build_archive.sh
#
# Collects every script used and every output produced across the whole
# project into one folder, with a manifest documenting what's there and what
# produced it. For the committee and for the permanent record.
#
# RUN (on Grace, from anywhere):
#   bash build_archive.sh
#
# OUTPUT:
#   /scratch/user/$USER/pdi/ARCHIVE_<date>/         the folder
#   /scratch/user/$USER/pdi_archive_<date>.tar.gz    the packaged version
#
# WHAT IT DOES NOT ASSUME
# -----------------------------------------------------------------------------
# Not every script that was written was necessarily run to completion, and not
# every run necessarily produced every possible output file (e.g. sPLS was
# blocked by a compiler issue on this cluster; a few sensitivity scripts were
# built but their run status was never explicitly confirmed). This script
# copies whatever it can FIND, does not fail if something is missing, and
# writes every miss into MANIFEST.md under "not found" rather than silently
# skipping it - so the manifest is an honest record of what exists, not an
# assumption that everything ran.
# =============================================================================
set -uo pipefail

ROOT="/scratch/user/$USER/pdi"
STAMP=$(date +%Y%m%d_%H%M)
ARCHIVE="${ROOT}/ARCHIVE_${STAMP}"
SCRIPTS_DIR="${ROOT}/scripts"
R_DIR="${SCRIPTS_DIR}/R_scripts"
GHI_DIR="${ROOT}/r_analysis/ghi"
GDI_DIR="${ROOT}/r_analysis/gdi"
GEN_DIR="${ROOT}/r_analysis/genus"

mkdir -p "${ARCHIVE}"/{01_scripts/bash,01_scripts/R,02_metadata,03_ghi_outputs,04_gdi_outputs,04_gdi_outputs/sensitivity_analyses,05_bioinformatics_outputs,06_figures,07_handover,08_qiime2_diagnostics}

FOUND_LOG="${ARCHIVE}/_found.txt"
MISSING_LOG="${ARCHIVE}/_missing.txt"
> "$FOUND_LOG"; > "$MISSING_LOG"

# copy_if src dest_dir description
copy_if() {
  local src="$1" dest="$2" desc="$3"
  if [ -e "$src" ]; then
    cp -rp "$src" "$dest"/ 2>/dev/null
    echo "${desc}|$(basename "$src")|$(du -h "$src" 2>/dev/null | cut -f1)" >> "$FOUND_LOG"
  else
    echo "${desc}|$(basename "$src")" >> "$MISSING_LOG"
  fi
}

echo "Building archive at ${ARCHIVE}"
echo ""

# -----------------------------------------------------------------------------
# 1. SCRIPTS - everything used, whether or not it produced output
# -----------------------------------------------------------------------------
echo "[1/8] Scripts"
for f in "${SCRIPTS_DIR}"/*.sh; do
  copy_if "$f" "${ARCHIVE}/01_scripts/bash" "QIIME2/bash stage"
done
for f in "${R_DIR}"/*.R; do
  copy_if "$f" "${ARCHIVE}/01_scripts/R" "R script"
done

# -----------------------------------------------------------------------------
# 2. METADATA - the corrected file, the pre-correction backups (audit trail),
#    and the QIIME2 sample metadata
# -----------------------------------------------------------------------------
echo "[2/8] Metadata (current + backups, as an audit trail of the classification corrections)"
copy_if "${SCRIPTS_DIR}/GDI_Metadata.xlsx"            "${ARCHIVE}/02_metadata" "current, corrected metadata"
for f in "${SCRIPTS_DIR}"/GDI_Metadata_backup_*.xlsx "${SCRIPTS_DIR}"/GDI_Metadata_CORRECTED.xlsx; do
  [ -e "$f" ] && copy_if "$f" "${ARCHIVE}/02_metadata" "pre-correction backup"
done
for f in "${SCRIPTS_DIR}"/PDI_metadata_qiime2*.tsv; do
  [ -e "$f" ] && copy_if "$f" "${ARCHIVE}/02_metadata" "QIIME2 sample metadata"
done
copy_if "${SCRIPTS_DIR}/treatment_aliases.tsv" "${ARCHIVE}/02_metadata" "superseded - aliases replaced by direct workbook correction"

# -----------------------------------------------------------------------------
# 3. GHI OUTPUTS
# -----------------------------------------------------------------------------
echo "[3/8] GHI outputs"
copy_if "${GHI_DIR}/ghi_scores_pca.tsv"          "${ARCHIVE}/03_ghi_outputs" "per-sample GHI scores - the index itself"
copy_if "${GHI_DIR}/ghi_index_definition.tsv"    "${ARCHIVE}/03_ghi_outputs" "portable index definition (for scoring new data)"
copy_if "${GHI_DIR}/ghi_pca_summary.txt"         "${ARCHIVE}/03_ghi_outputs" "domain construction log (04b)"
copy_if "${GHI_DIR}/mixed_model_validation.txt"  "${ARCHIVE}/03_ghi_outputs" "primary GHI validation (04e) - additive/challenge contrasts"
copy_if "${GHI_DIR}/ghi_convergent_validity.txt" "${ARCHIVE}/03_ghi_outputs" "morphology vs biomarker correlation check (04f)"
copy_if "${GHI_DIR}/ghi_domain_correlation.txt"  "${ARCHIVE}/03_ghi_outputs" "domain correlation matrix (04g)"
copy_if "${GHI_DIR}/ghi_ordinal.txt"             "${ARCHIVE}/03_ghi_outputs" "single-scale challenge/control/additive test (15)"
copy_if "${GHI_DIR}/standardisation_comparison.txt" "${ARCHIVE}/03_ghi_outputs" "z-scoring sensitivity check"
copy_if "${GHI_DIR}/validation_matrix.txt"       "${ARCHIVE}/03_ghi_outputs" "superseded - concordance analysis, replaced by 04e"

# -----------------------------------------------------------------------------
# 4. GDI OUTPUTS - primary
# -----------------------------------------------------------------------------
echo "[4/8] GDI outputs - primary"
copy_if "${GDI_DIR}/gdi_panel.txt"              "${ARCHIVE}/04_gdi_outputs" "final panel (empty at full-corpus, 70% threshold)"
copy_if "${GDI_DIR}/gdi_stability.tsv"          "${ARCHIVE}/04_gdi_outputs" "per-genus stability selection frequencies"
copy_if "${GDI_DIR}/gdi_index_definition.tsv"   "${ARCHIVE}/04_gdi_outputs" "portable GDI definition"
copy_if "${GDI_DIR}/gdi_scores.tsv"             "${ARCHIVE}/04_gdi_outputs" "per-arm GDI scores"
copy_if "${GDI_DIR}/gdi_construct_summary.txt"  "${ARCHIVE}/04_gdi_outputs" "GDI construction log (07)"
copy_if "${GDI_DIR}/lasso_summary.txt"          "${ARCHIVE}/04_gdi_outputs" "LASSO fit summary (06)"
copy_if "${GDI_DIR}/loso_summary.txt"           "${ARCHIVE}/04_gdi_outputs" "leave-one-study-out validation (08)"
copy_if "${GDI_DIR}/loso_predictions.tsv"       "${ARCHIVE}/04_gdi_outputs" "per-arm held-out predictions"
copy_if "${GDI_DIR}/loso_per_study.tsv"         "${ARCHIVE}/04_gdi_outputs" "per-study LOSO correlations"
copy_if "${GDI_DIR}/maaslin2_all_results.tsv"   "${ARCHIVE}/04_gdi_outputs" "MaAsLin2 screening, all 84 genera"
copy_if "${GDI_DIR}/maaslin2_selected_genera.txt" "${ARCHIVE}/04_gdi_outputs" "genera passing FDR screen"
copy_if "${GDI_DIR}/maaslin2_summary.txt"       "${ARCHIVE}/04_gdi_outputs" "MaAsLin2 run log"
copy_if "${GDI_DIR}/absolute_discriminant_test.txt" "${ARCHIVE}/04_gdi_outputs" "single-sample deployability test (09)"
copy_if "${GDI_DIR}/variance_partition.txt"     "${ARCHIVE}/04_gdi_outputs" "PERMANOVA: age vs treatment on community structure"
copy_if "${GDI_DIR}/ghi_permanova.txt"          "${ARCHIVE}/04_gdi_outputs" "PERMANOVA: does GHI explain community structure (13)"
copy_if "${GDI_DIR}/taxon_metaanalysis.txt"     "${ARCHIVE}/04_gdi_outputs" "per-taxon random-effects meta-analysis (18)"
copy_if "${GDI_DIR}/taxon_metaanalysis.tsv"     "${ARCHIVE}/04_gdi_outputs" "per-taxon meta-analysis, machine-readable"

# with-covariates MaAsLin2 runs (naming may vary by which flags were set)
for d in "${GDI_DIR}"/maaslin2_age* "${GDI_DIR}"/maaslin2; do
  [ -d "$d" ] && copy_if "$d" "${ARCHIVE}/04_gdi_outputs" "MaAsLin2 output directory (check contents for which covariates)"
done

# -----------------------------------------------------------------------------
# 5. GDI OUTPUTS - sensitivity analyses and alternative constructions
# -----------------------------------------------------------------------------
echo "[5/8] GDI sensitivity analyses"
S="${ARCHIVE}/04_gdi_outputs/sensitivity_analyses"
copy_if "${GDI_DIR}/lambda_comparison"          "$S" "lambda.1se vs lambda.min (06c)"
copy_if "${GDI_DIR}/stratified_test.txt"        "$S" "all/additive/challenge/age-adjusted arm subsets (06e)"
copy_if "${GDI_DIR}/challenge_gdi.txt"          "$S" "challenge-only supervised panel + LOSO (06f)"
copy_if "${GDI_DIR}/challenge_gdi_predictions.tsv" "$S" "challenge-only LOSO predictions"
copy_if "${GDI_DIR}/gdi_by_class.txt"           "$S" "GDI fitted separately within each class (14)"
copy_if "${GDI_DIR}/unsupervised_gdi.txt"       "$S" "unsupervised indices, control-based (10)"
copy_if "${GDI_DIR}/unsupervised_scores.tsv"    "$S" "unsupervised index scores"
copy_if "${GDI_DIR}/axis_loso.txt"              "$S" "control-challenge axis, leave-one-study-out (11)"
copy_if "${GDI_DIR}/axis_loso_scores.tsv"       "$S" "axis LOSO scores"
copy_if "${GDI_DIR}/gdi_three_ways.txt"         "$S" "deviation/axis/fixed-panel comparison (12)"
copy_if "${GDI_DIR}/gdi_three_ways_scores.tsv"  "$S" "three-ways per-arm scores"
copy_if "${GDI_DIR}/three_way_axis.txt"         "$S" "control-challenge-additive axis (16)"
copy_if "${GDI_DIR}/three_way_axis_scores.tsv"  "$S" "three-way axis scores"
copy_if "${GDI_DIR}/panel_size_test.txt"        "$S" "panel size N=3..10 comparison (06d), if run"
copy_if "${GDI_DIR}/spls"                       "$S" "sparse PLS output (19), if the cluster compiler issue was resolved"
copy_if "${GDI_DIR}/maaslin2_class"             "$S" "MaAsLin2 - challenge/additive vs control (17)"

# -----------------------------------------------------------------------------
# 6. BIOINFORMATICS / GENUS-LEVEL OUTPUTS
# -----------------------------------------------------------------------------
echo "[6/8] Genus-level and batch-correction outputs"
copy_if "${GEN_DIR}/genus_counts.tsv"                  "${ARCHIVE}/05_bioinformatics_outputs" "raw genus counts, post-prevalence-filter"
copy_if "${GEN_DIR}/genus_proportions_batch_adjusted.tsv" "${ARCHIVE}/05_bioinformatics_outputs" "MMUPHin-corrected proportions"
copy_if "${GEN_DIR}/genus_clr_batch_adjusted.tsv"      "${ARCHIVE}/05_bioinformatics_outputs" "final CLR matrix used by the GDI"
copy_if "${GEN_DIR}/mmuphin_diagnostic.pdf"            "${ARCHIVE}/05_bioinformatics_outputs" "MMUPHin's own shrinkage/adjustment diagnostic"
copy_if "${GEN_DIR}/cmultRepl_dropped_samples.txt"     "${ARCHIVE}/05_bioinformatics_outputs" "samples dropped at zero-replacement (>80% zeros)"
copy_if "${GEN_DIR}/cmultRepl_dropped_genera.txt"      "${ARCHIVE}/05_bioinformatics_outputs" "genera dropped at zero-replacement"

# -----------------------------------------------------------------------------
# 7. FIGURES
# -----------------------------------------------------------------------------
echo "[7/8] Figures"
copy_if "${ROOT}/figures" "${ARCHIVE}/06_figures" "deck-ready figures (make_figures.R)"

# -----------------------------------------------------------------------------
# 8. HANDOVER PACKAGE + QIIME2-STAGE DIAGNOSTICS
# -----------------------------------------------------------------------------
echo "[8/8] Handover package and QIIME2 diagnostics"
copy_if "${ROOT}/handover" "${ARCHIVE}/07_handover" "standalone scoring package (index definitions + scripts + examples)"
copy_if "${ROOT}/merged/braycurtis_emperor.qzv"                    "${ARCHIVE}/08_qiime2_diagnostics" "PCoA (all 50 studies, pre-restriction)"
copy_if "${ROOT}/merged/jaccard_emperor.qzv"                       "${ARCHIVE}/08_qiime2_diagnostics" "PCoA, presence/absence"
copy_if "${ROOT}/merged/braycurtis_Study_ID_permanova.qzv"         "${ARCHIVE}/08_qiime2_diagnostics" "PERMANOVA - batch effect"
copy_if "${ROOT}/merged/braycurtis_Intervention_Class_permanova.qzv" "${ARCHIVE}/08_qiime2_diagnostics" "PERMANOVA - biological signal"

# =============================================================================
# MANIFEST
# =============================================================================
N_FOUND=$(wc -l < "$FOUND_LOG")
N_MISSING=$(wc -l < "$MISSING_LOG")

cat > "${ARCHIVE}/MANIFEST.md" << MANIFEST_EOF
# Project Archive
Built $(date "+%Y-%m-%d %H:%M")
PROSPERO CRD420261384137

## What this is
Every script used and every output produced across the full pipeline -
systematic review corpus through GHI/GDI construction, validation, and every
sensitivity analysis run. Organized so a committee member or a future reader
can trace any number in the manuscript back to the file and script that
produced it.

## Folder structure

- \`01_scripts/\`       every bash (QIIME2) and R script used, in the order they run
- \`02_metadata/\`      the corrected metadata workbook, plus every backup made
                       before a correction was applied - this is the audit
                       trail for the intervention-classification corrections
- \`03_ghi_outputs/\`   Gut Health Index construction and validation
- \`04_gdi_outputs/\`   Gut Dysbiosis Index construction and validation, split
                       into primary results and \`sensitivity_analyses/\` -
                       every arm subset, target, and construction tested
- \`05_bioinformatics_outputs/\`  genus-level tables and batch-correction diagnostics
- \`06_figures/\`       deck-ready figures
- \`07_handover/\`      standalone package for applying either index to new data
- \`08_qiime2_diagnostics/\`  PCoA and PERMANOVA from the QIIME2 stage

## Files found: ${N_FOUND}
## Files/paths not found: ${N_MISSING}
(see _found.txt and _missing.txt for the itemised lists - a "not found" entry
usually means that particular sensitivity script was written but its run was
never confirmed to completion, or produces output only under a flag that
wasn't set this session. It does not necessarily mean the analysis failed.)

## Reading order for the numbers that matter most

1. \`03_ghi_outputs/mixed_model_validation.txt\` - the primary GHI validation
2. \`03_ghi_outputs/ghi_ordinal.txt\` - the single-scale confirmatory test
3. \`04_gdi_outputs/gdi_stability.tsv\` - why no full-corpus panel exists
4. \`04_gdi_outputs/ghi_permanova.txt\` - the direct test of whether GHI
   explains any community variance (it does not, at full corpus)
5. \`04_gdi_outputs/sensitivity_analyses/challenge_gdi.txt\` - the
   challenge-only subset, the closest anything came to a working panel, and
   why it did not survive honest cross-validation
6. \`04_gdi_outputs/taxon_metaanalysis.txt\` - which genera are consistent
   across studies independent of the LASSO panel

## Known caveats to carry into the report

- The classification correction (see \`02_metadata/\` backups) changed which
  arms count as challenge vs additive; every GDI sensitivity analysis in
  \`sensitivity_analyses/\` should be read as POST-correction unless its
  filename/timestamp says otherwise.
- \`spls/\` and \`panel_size_test.txt\` may be absent - sPLS was blocked by a
  compiler mismatch on this cluster and the panel-size test's run was never
  explicitly confirmed.
MANIFEST_EOF

echo ""
echo "=== DONE ==="
echo "Archive: ${ARCHIVE}"
echo "Files found:   ${N_FOUND}"
echo "Files missing: ${N_MISSING}   (see ${MISSING_LOG})"
echo ""
echo "To package:"
echo "  cd ${ROOT}"
echo "  tar -czf ${ROOT}/../pdi_archive_${STAMP}.tar.gz -C ${ROOT} $(basename "$ARCHIVE")"
echo "  (or: tar -czf ~/pdi_archive_${STAMP}.tar.gz -C ${ROOT} $(basename "$ARCHIVE"))"
