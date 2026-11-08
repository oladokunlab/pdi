#!/usr/bin/env bash
# =============================================================================
# organize_final.sh (v2 - after the full corrected rebuild)
#
# Four folders:
#   A_Pipeline_Scripts/      only the scripts that ARE the guide's 25 steps
#   B_Generated_Files/       the outputs those scripts produced
#   C_Validation_and_Checks/ scripts/outputs necessary because an issue was
#                            found, but not steps in the guide
#   D_QC_Reports/            placeholder - FastQC/MultiQC added manually
#
# Superseded by this version: gdi_scores.tsv is NOT copied. The final,
# fully-corrected pipeline produces only 2 genera clearing the panel floor
# (registered minimum is 3), so 07_gdi_construct.R correctly refuses to build
# a centroid. There is no full-corpus panel to archive - this is the result,
# not a missing file. See README.
#
# RUN: bash organize_final.sh
# =============================================================================
set -uo pipefail

ROOT="/scratch/user/$USER/pdi"
STAMP=$(date +%Y%m%d_%H%M)
OUT="${ROOT}/FINAL_${STAMP}"
SCR="${ROOT}/scripts"
R="${SCR}/R_scripts"
GHI="${ROOT}/r_analysis/ghi"
GDI="${ROOT}/r_analysis/gdi"
GEN="${ROOT}/r_analysis/genus"

mkdir -p "${OUT}"/A_Pipeline_Scripts/{bash,R} \
         "${OUT}"/B_Generated_Files/{GHI,GDI,genus_level} \
         "${OUT}"/C_Validation_and_Checks/{scripts,outputs} \
         "${OUT}"/D_QC_Reports \
         "${OUT}"/E_Figures

cp_if() { [ -e "$1" ] && cp -rp "$1" "$2"/ 2>/dev/null; }

# =============================================================================
# A. PIPELINE SCRIPTS
# =============================================================================
A="${OUT}/A_Pipeline_Scripts"

cp_if "${SCR}/00_config.sh"                 "$A/bash"
cp_if "${SCR}/01_build_manifests.sh"        "$A/bash"
cp_if "${SCR}/02_denoise_per_study.sh"      "$A/bash"
cp_if "${SCR}/02b_denoise_summary.sh"       "$A/bash"
cp_if "${SCR}/03_merge_and_taxonomy.sh"     "$A/bash"
cp_if "${SCR}/04_filter_and_export.sh"      "$A/bash"

cp_if "${R}/00_r_config.R"                  "$A/R"
cp_if "${R}/01_genus_agglomeration.R"       "$A/R"
cp_if "${R}/03_batch_correction.R"          "$A/R"
cp_if "${R}/04b_ghi_pca.R"                  "$A/R"      # Performance split, morphology fix
cp_if "${R}/04e_mixed_model_validation.R"   "$A/R"      # domain list synced with 04b
cp_if "${R}/05_maaslin2.R"                  "$A/R"      # AGE_COVARIATE + CLASS_COVARIATE + SEGMENT_COVARIATE
cp_if "${R}/19_sparse_pls.R"                "$A/R"      # attempted (mixOmics) - see C for the working fallback
cp_if "${R}/18_taxon_metaanalysis.R"        "$A/R"
cp_if "${R}/12_gdi_three_ways.R"            "$A/R"
cp_if "${R}/06_lasso_selection.R"           "$A/R"
cp_if "${R}/07_gdi_construct.R"             "$A/R"      # ran; produced no panel, see README
cp_if "${R}/08_loso_validate.R"             "$A/R"
cp_if "${R}/04g_domain_correlation.R"       "$A/R"      # generates the domain heatmap input
cp_if "${R}/make_figures.R"                 "$A/R"

# =============================================================================
# B. GENERATED FILES
# =============================================================================
B="${OUT}/B_Generated_Files"

cp_if "${GHI}/ghi_scores_pca.tsv"            "$B/GHI"
cp_if "${GHI}/ghi_index_definition.tsv"      "$B/GHI"
cp_if "${GHI}/mixed_model_validation.txt"    "$B/GHI"
cp_if "${GHI}/ghi_pca_summary.txt"           "$B/GHI"
cp_if "${GHI}/ghi_ordinal.txt"               "$B/GHI"
cp_if "${GHI}/ghi_domain_correlation.txt"        "$B/GHI"
cp_if "${GHI}/ghi_domain_correlation_matrix.tsv" "$B/GHI"
cp_if "${GHI}/ghi_diagnostics.txt"               "$B/GHI"
cp_if "${GHI}/microbiome_as_domain.txt"          "$B/GHI"
cp_if "${GHI}/domain_restructure.txt"            "$B/GHI"

cp_if "${GDI}/gdi_stability.tsv"             "$B/GDI"
cp_if "${GDI}/maaslin2_all_results.tsv"      "$B/GDI"
cp_if "${GDI}/maaslin2_selected_genera.txt"  "$B/GDI"
cp_if "${GDI}/taxon_metaanalysis.txt"        "$B/GDI"
cp_if "${GDI}/taxon_metaanalysis.tsv"        "$B/GDI"
cp_if "${GDI}/loso_summary.txt"              "$B/GDI"
cp_if "${GDI}/loso_predictions.tsv"          "$B/GDI"
cp_if "${GDI}/maaslin2_age_class_seg"        "$B/GDI"   # final covariate-extended MaAsLin2 run
cp_if "${GDI}/gdi_model.rds"                 "$B/GDI"   # was missing from v1 - the underlying LASSO/stability object
# gdi_scores.tsv deliberately NOT copied - no panel was constructed on the
# final corrected data (2 genera cleared the floor, 3 required). See README.

cp_if "${GEN}/genus_counts.tsv"                     "$B/genus_level"
cp_if "${GEN}/genus_clr_batch_adjusted.tsv"         "$B/genus_level"
cp_if "${GEN}/genus_proportions_batch_adjusted.tsv" "$B/genus_level"
cp_if "${GEN}/batch_correction_summary.txt"         "$B/genus_level"
cp_if "${GEN}/genus_agglomeration_summary.txt"      "$B/genus_level"
cp_if "${GEN}/cmultRepl_dropped_genera.txt"         "$B/genus_level"   # was missing from v1
cp_if "${GEN}/cmultRepl_dropped_samples.txt"        "$B/genus_level"   # was missing from v1

# =============================================================================
# C. VALIDATION AND CHECKS
# =============================================================================
C="${OUT}/C_Validation_and_Checks"

cp_if "${R}/fix_classes.R"                  "$C/scripts"
cp_if "${R}/fix_classes2.R"                 "$C/scripts"
cp_if "${R}/morph_check.R"                  "$C/scripts"
cp_if "${R}/units_audit.R"                  "$C/scripts"
cp_if "${R}/barrier_and_coverage.R"         "$C/scripts"
cp_if "${R}/variance_partition.R"           "$C/scripts"
cp_if "${R}/13_ghi_permanova.R"             "$C/scripts"
cp_if "${R}/02b_batch_confounding_check.R"  "$C/scripts"   # diagnosed the Segment/Study_ID confound
cp_if "${R}/06e_stratified_test.R"          "$C/scripts"
cp_if "${R}/06f_challenge_gdi.R"            "$C/scripts"
cp_if "${R}/14_gdi_by_class.R"              "$C/scripts"
cp_if "${R}/16_three_way_axis.R"            "$C/scripts"
cp_if "${R}/17_maaslin_class.R"             "$C/scripts"
cp_if "${R}/20_control_vs_treated.R"        "$C/scripts"
cp_if "${R}/21_target_comparison_domains.R" "$C/scripts"
cp_if "${R}/19b_sparse_pls_fallback.R"      "$C/scripts"   # the sPLS method that actually ran
cp_if "${R}/22_binary_presence_lasso.R"     "$C/scripts"   # alternative method: binary presence/absence
cp_if "${R}/23_percentile_normalization.R"  "$C/scripts"   # alternative method: percentile normalisation
cp_if "${R}/24_elastic_net.R"               "$C/scripts"   # alternative method: elastic net alpha sweep
cp_if "${R}/25_random_forest.R"             "$C/scripts"   # alternative method: random forest + LOSO
cp_if "${R}/fix_missing_ages.R"             "$C/scripts"   # age audit trail
cp_if "${R}/fix_ages_final.R"               "$C/scripts"   # age audit trail
cp_if "${R}/fix_ma2024_split.R"             "$C/scripts"   # Ma_2024 two-experiment split
cp_if "${R}/rename_ma2024_back.R"           "$C/scripts"   # Ma_2024 naming
cp_if "${R}/fix_arm_names.R"                "$C/scripts"   # arm-name standardisation
cp_if "${R}/verify_metadata_full.R"         "$C/scripts"   # reusable validator
cp_if "${R}/verify_ages_vs_sources.R"    "$C/scripts"   # age verification against sources
cp_if "${R}/26_microbiome_as_domain.R"      "$C/scripts"   # microbiome tested as a GHI domain
cp_if "${R}/27_runlevel_gdi.R"              "$C/scripts"   # run-level classifier, treatment label
cp_if "${R}/28_runlevel_groupings.R"        "$C/scripts"   # alternative treatment groupings
cp_if "${R}/29_ghi_diagnostics.R"           "$C/scripts"   # six GHI diagnostics
cp_if "${R}/30_axis_vs_highcoverage_ghi.R"  "$C/scripts"   # axis vs coverage-restricted GHI
cp_if "${R}/31_picrust_scfa_check.R"        "$C/scripts"   # PICRUSt2 validated against measured SCFA
cp_if "${R}/32_bactotraits.R"               "$C/scripts"   # trait CWM and Rao's Q
cp_if "${R}/33_guild_ratios.R"              "$C/scripts"   # functional guild ratios
cp_if "${R}/34_domain_restructure.R"        "$C/scripts"   # domain merge tests
cp_if "${R}/35_scfa_diagnosis.R"            "$C/scripts"   # why the SCFA domain inverts
cp_if "${R}/36_siga_orientation.R"          "$C/scripts"   # sIgA direction and Ig merge test
cp_if "${R}/audit_compartments.R"           "$C/scripts"   # serum vs mucosa audit
cp_if "${R}/check_ig_compartment.R"         "$C/scripts"   # immunoglobulin compartment check
cp_if "${R}/fix_iga_compartment.R"          "$C/scripts"   # 64 serum IgA reassigned
cp_if "${R}/fix_siga_mislabel.R"            "$C/scripts"   # 15 mislabelled sIgA reassigned
cp_if "${R}/fix_zhang2026a_revert.R"        "$C/scripts"   # reverted a wrong age correction
cp_if "${R}/diagnose_bacto.R"               "$C/scripts"   # trait matrix rank diagnosis
cp_if "${SCR}/05_picrust2.sh"               "$C/scripts"   # PICRUSt2 pipeline
cp_if "${SCR}/../verify_dependencies.sh"    "$C/scripts" 2>/dev/null
cp_if "${SCR}/verify_dependencies.sh"       "$C/scripts"

cp_if "${GDI}/variance_partition.txt"                 "$C/outputs"
cp_if "${GDI}/ghi_permanova.txt"                       "$C/outputs"
cp_if "${GDI}/stratified_test.txt"                     "$C/outputs"
cp_if "${GDI}/challenge_gdi.txt"                       "$C/outputs"
cp_if "${GDI}/gdi_by_class.txt"                        "$C/outputs"
cp_if "${GDI}/three_way_axis.txt"                      "$C/outputs"
cp_if "${GDI}/control_vs_treated.txt"                  "$C/outputs"
cp_if "${GDI}/gdi_three_ways.txt"                      "$C/outputs"   # was missing - output of 12, script was archived but not its result
cp_if "${GDI}/gdi_three_ways_scores.tsv"                "$C/outputs"  # same
cp_if "${GDI}/target_comparison_domains.txt"           "$C/outputs"
cp_if "${GDI}/binary_presence_lasso.txt"               "$C/outputs"
cp_if "${GDI}/percentile_normalization.txt"            "$C/outputs"
cp_if "${GDI}/elastic_net_sweep.txt"                   "$C/outputs"
cp_if "${GDI}/random_forest.txt"                       "$C/outputs"
cp_if "${GDI}/runlevel_gdi.txt"                        "$C/outputs"
cp_if "${GDI}/runlevel_groupings.txt"                  "$C/outputs"
cp_if "${GDI}/bactotraits.txt"                         "$C/outputs"
cp_if "${GDI}/guild_ratios.txt"                        "$C/outputs"
cp_if "${GDI}/axis_vs_coverage.txt"                    "$C/outputs"
cp_if "${GDI}/spls_fallback"                           "$C/outputs"
cp_if "${GDI}/maaslin2_class"                          "$C/outputs"
cp_if "${SCR}/PDI_metadata_qiime2_FINAL_v5.tsv"          "$C/outputs" 2>/dev/null
cp_if "${SCR}/GDI_metadata_qiime2.tsv"                   "$C/outputs" 2>/dev/null
cp_if "${SCR}/GDI_Metadata_backup_20260818_005055.xlsx" "$C/outputs"
cp_if "${SCR}/GDI_Metadata_backup_20260819_152409.xlsx" "$C/outputs"
cp_if "${SCR}/GDI_Metadata.xlsx"                        "$C/outputs"

# =============================================================================
# G. PICRUSt2 - functional prediction
# Only the summary tables are archived. The full PICRUSt2 output directory is
# several GB and regenerable from 05_picrust2.sh plus the exported ASV
# sequences, so it is left in place rather than copied.
# =============================================================================
mkdir -p "${OUT}/G_PICRUSt2"
cp_if "${ROOT}/picrust2/out/KO_metagenome_out/pred_metagenome_unstrat.tsv" "${OUT}/G_PICRUSt2"
cp_if "${ROOT}/picrust2/out/pathways_out/path_abun_unstrat.tsv"            "${OUT}/G_PICRUSt2"
cp_if "${ROOT}/picrust2/out/KO_metagenome_out/weighted_nsti.tsv"           "${OUT}/G_PICRUSt2"

# =============================================================================
# F. QIIME2 DIAGNOSTICS - refreshed after the metadata corrections
# =============================================================================
mkdir -p "${OUT}/F_QIIME2_Diagnostics"
for f in braycurtis_emperor.qzv jaccard_emperor.qzv merged-table-summary.qzv \
         braycurtis_Study_ID_permanova.qzv \
         braycurtis_Intervention_Class_permanova.qzv \
         braycurtis_PCoA_Group_permanova.qzv; do
  cp_if "${ROOT}/merged/${f}" "${OUT}/F_QIIME2_Diagnostics"
done
cp_if "${ROOT}/merged/pcoa_depth_used.txt" "${OUT}/F_QIIME2_Diagnostics"
cp_if "${SCR}/03a_merge_and_pcoa.sh"           "$A/bash"
cp_if "${SCR}/03b_build_phylogeny.sh"          "$A/bash"
cp_if "${SCR}/03c_refresh_pcoa_permanova.sh"   "$A/bash"

# =============================================================================
# E. FIGURES
# =============================================================================
cp_if "${ROOT}/figures" "${OUT}/E_Figures"

# =============================================================================
# D. QC REPORTS - placeholder, populated manually
# =============================================================================
cat > "${OUT}/D_QC_Reports/PUT_MULTIQC_REPORT_HERE.txt" << 'EOF'
Add here:
  - multiqc_report.html (pooled, all 50 studies)
  - per_study/<Study>/<Study>_multiqc.html  (if generated)

Source: /scratch/user/$USER/pdi/qc_reports/
EOF

# =============================================================================
# README
# =============================================================================
cat > "${OUT}/README.md" << 'MDEOF'
# PDI Project - Final Organized Set (v2, post full-corpus rebuild)

Rebuilt after: the Performance/Biomarker split, correcting a real MMUPHin
error (age was confounded with study identity in 30 of 48 studies and had to
be reverted there, while remaining a valid MaAsLin2 covariate), and a full,
verified rebuild of the entire pipeline. Every file here was checked against
its inputs with verify_dependencies.sh before this archive was built.

## Corpus, as of this rebuild
115 genera x 2064 samples across 44 V3-V4 studies survive prevalence
filtering (post CLR zero-replacement counts are in genus_level/
batch_correction_summary.txt). This differs from the 84-genera figure seen
earlier in the project; the increase is a direct, verified consequence of the
corrected classification, not an error - traced and confirmed before this
archive was built.

## A_Pipeline_Scripts/
The 25-step pipeline. One methodological correction made during this
rebuild, stated precisely:

  Age (Age_d) was tested as a protected covariate in MMUPHin's adjust_batch
  and found to be mathematically confounded with study identity: 30 of 48
  studies sample birds at a single age, making "this study" and "this age"
  the same information for adjust_batch's purposes, and it correctly refuses
  to run rather than produce an ill-defined adjustment. Age was REMOVED from
  MMUPHin (03_batch_correction.R now protects Intervention_Class only) and
  RETAINED in MaAsLin2 (05_maaslin2.R, where study enters as a random effect
  and this constraint does not apply). This is why 03 and 05 are run with
  different covariate flags - see the run commands recorded in
  C_Validation_and_Checks/scripts/verify_dependencies.sh's companion notes.

Two further, already-documented divergences from the guide (LASSO centering,
study-level bootstrap) are unchanged from earlier and remain in 06.

## B_Generated_Files/
Direct output of the scripts in A.

IMPORTANT: gdi_scores.tsv is NOT included. On the fully corrected pipeline,
07_gdi_construct.R found only 2 genera clearing the 20% stability floor
(3 is the minimum required to build a centroid) and correctly halted rather
than force a panel. This is the result, not a missing file - see
gdi_stability.tsv for the 2 genera and their stability (Peptococcus 26%,
Limosilactobacillus 24%).

## C_Validation_and_Checks/
Includes the full audit trail: both classification-correction scripts with
their metadata backups, the morphology-bug diagnostic, the MMUPHin
confounding diagnostic (02b), and every sensitivity analysis - arm subsets,
per-class fits, unsupervised axes (signed and pooled), the GHI-free
class-level MaAsLin2 test, the four-target domain comparison, and the
sparse PLS result (spls package - the mixOmics route failed after three
independent infrastructure problems on this cluster: a compiler flag
mismatch, a session-lingering issue, and a Lustre filesystem issue unpacking
Boost headers; 19_sparse_pls.R is kept as the attempted, unrun script for
the record, 19b_sparse_pls_fallback.R is the one that actually produced a
result).

## D_QC_Reports/
Placeholder. FastQC/MultiQC reports to be added manually from
/scratch/user/$USER/pdi/qc_reports/.

## E_Figures/
Deck figures from make_figures.R. slide23_fold_consistency.png will be
absent - correctly - since no panel exists to show fold consistency for.
slide16_domain_heatmap.png now shows 10 domains (Performance, Mucin and
Anti_Inflammatory added since the original 7-domain figure was designed).

## F_QIIME2_Diagnostics/
PCoA plots and PERMANOVA tests, regenerated after the metadata corrections
(03c_refresh_pcoa_permanova.sh). The distance matrices, ordination
coordinates and phylogeny were NOT recomputed - they derive from sequences
and counts, which the metadata corrections did not touch. Only the
metadata-dependent visualisations and group tests were rebuilt. The
Study_ID PERMANOVA in particular differs from the pre-correction version,
since Ma_2024 is now a single correctly-scoped study rather than two
conflated experiments.

## Metadata corrections applied in this round
- 94 Age_d values corrected across 11 studies, each verified against the
  source publication (Zhao_2023a tissue timepoint, Massacci_2019 challenge
  day mistaken for a sampling day, Zhao_2022 performance-only endpoint,
  Zhang_2026a broiler defaults in a 180-day layer study, Jiang_2023
  days-post-infection stored as absolute age, and six smaller corrections).
- Ma_2024 separated into its two underlying experiments; only the FMT chick
  trial has sequencing data, and its age corrected from 21 to 18 days.
- 190 sequencing runs across four studies relabelled so treatment-arm names
  match between the 16S metadata and the biomarker workbook.

## Known limitations, disclosed rather than fixed
- Jansseune_2024: the public SRA deposit records diet type only, never the
  additive condition (Ctrl/Pro/Post). Confirmed against the run table and
  the paper. Its 162 runs are handled by study-level GHI averaging.
- Zhao_2022: three rows have no intestinal segment recorded in the source.

## The index is SEVEN domains, not nine
Morphology, Barrier (ZO-1, occludin, claudin-1, MUC2), Inflammation,
Anti_Inflammatory, SCFA, Performance, Immunoglobulin_Systemic (IgG, IgM,
serum IgA). 18 markers, 529 samples.

Two structural changes were made after the initial build, both on biology
first and checked against the corpus afterwards:

  MUC2 folded into Barrier. Mucus is part of the intestinal barrier and the
  poultry literature reports it that way; it correlated +0.754 with Barrier,
  higher than Barrier with Morphology. As a one-marker domain it also carried
  the same weight in the average as Morphology's 311 samples.

  Immunoglobulin_Mucosal (sIgA) removed from the SCORED index and reported
  separately. With its compartment cleaned it scores AUC 0.071 - sharply
  inverted - but the sign cannot be flipped, because sIgA rises under
  challenge AND under beneficial additives. It marks mucosal immune activity,
  not the direction of health. SCORE_SIGA=yes in 04b restores it.

Two candidate merges were tested and rejected: inflammation with
anti-inflammatory (they rise together in real data, so averaging cancels real
signal) and the two immunoglobulin domains (circulating Ig marks the insult
rather than health).

## Compartment corrections
An audit of every dual-compartment biomarker found two errors:

  All 64 IgA values in the corpus were measured in blood serum, confirmed
  against each source paper, while sitting in a domain labelled mucosal.
  Moved to Immunoglobulin_Systemic as IgA_serum.

  15 further values in the sIgA column came from three studies that report no
  secretory IgA at all (Jansseune_2024, Han_2023, Chen_2024). Also serum.

55 genuinely mucosal sIgA rows remain across 8 studies, all jejunal or ileal
homogenate. See audit_compartments.R and the two fix scripts in
C_Validation_and_Checks/scripts.

## One correction was wrong and was reverted
Zhang_2026a was changed to 180 days on the basis of a source quotation that
belonged to a DIFFERENT paper - the literature tool's document numbering did
not match this project's study IDs. Caught by checking the ENA study title for
the accession in this project's own metadata (PRJNA1223836), which states the
sampling ages outright: 21 and 42 days. All 41 ENA-resolvable studies were
subsequently reconciled against their archive titles with no further
mismatches. A citation is not evidence that the right document was read.

## Alternative methods tested after the registered method produced no panel
Binary presence/absence encoding, percentile normalisation against
within-study controls, an elastic-net alpha sweep, and random forest with
leave-one-study-out. All are in C_Validation_and_Checks. Every one returned
R2 at or near zero; these are reported findings, not failed attempts.

## Fixed in this rebuild (v2 gaps found during verification)
gdi_model.rds, cmultRepl_dropped_genera.txt and cmultRepl_dropped_samples.txt
were missing from the previous archive pass. All three are included now.

## Reading order for the numbers that matter most
1. B_Generated_Files/GHI/mixed_model_validation.txt - primary GHI validation
2. B_Generated_Files/GHI/ghi_ordinal.txt - confirmatory single-scale test
3. B_Generated_Files/GDI/gdi_stability.tsv - why no full-corpus panel exists
4. C_Validation_and_Checks/outputs/ghi_permanova.txt - direct test, no model built
5. C_Validation_and_Checks/outputs/challenge_gdi.txt - the classification-
   contamination catch (81%/AUC 0.638 before correction -> 61%/AUC 0.510 after)
6. C_Validation_and_Checks/outputs/spls_fallback/ - the Blautia/Enterococcus
   multivariate axis, the most coherent single finding in the project
   (n=31, hypothesis-generating)
MDEOF

N_A=$(find "${OUT}/A_Pipeline_Scripts" -type f | wc -l)
N_B=$(find "${OUT}/B_Generated_Files" -type f | wc -l)
N_C=$(find "${OUT}/C_Validation_and_Checks" -type f | wc -l)
N_E=$(find "${OUT}/E_Figures" -type f | wc -l)
N_F=$(find "${OUT}/F_QIIME2_Diagnostics" -type f 2>/dev/null | wc -l)
N_G=$(find "${OUT}/G_PICRUSt2" -type f 2>/dev/null | wc -l)

echo "=== DONE ==="
echo "A_Pipeline_Scripts:      ${N_A} files"
echo "B_Generated_Files:       ${N_B} files"
echo "C_Validation_and_Checks: ${N_C} files"
echo "D_QC_Reports:            placeholder only - add FastQC/MultiQC yourself"
echo "E_Figures:               ${N_E} files"
echo "F_QIIME2_Diagnostics:    ${N_F} files"
echo "G_PICRUSt2:              ${N_G} files"
echo ""
echo "Folder: ${OUT}"
echo ""
echo "To package:"
echo "  tar -czf ~/pdi_final_${STAMP}.tar.gz -C ${ROOT} $(basename "$OUT")"
