#!/bin/bash
# Rewrites r_analysis/gdi/maaslin2_summary.txt with the corrected text.
set -euo pipefail
D=/scratch/user/$USER/pdi/r_analysis/gdi
[ -f "$D/maaslin2_summary.txt" ] && mv "$D/maaslin2_summary.txt" "$D/maaslin2_summary_SUPERSEDED.txt"
cat > "$D/maaslin2_summary.txt" <<'PDIEOF'
=== MaAsLin2: genus associations with the GHI ===

Target (pre-specified): GHI_Combined
Runs analysed: 1625
Genera tested: 84
Studies contributing: 42  (43 studies are present in the CLR matrix; one contributes
no run that joins a GHI value, so 42 enter the model)

Pass FDR <= 0.25 : 28  (proceed to LASSO)
Pass FDR <= 0.05 : 10

Model AS RUN: per-genus linear mixed model with GHI as the fixed effect of
interest, ADJUSTED for Age_d, Challenge_Type and Segment as additional fixed
effects, and study identity as a random effect. Input already CLR-transformed,
so MaAsLin2 normalisation and transformation are both set to NONE.

  *** CORRECTION (see note 4 below). An earlier version of this file described
  *** the model as "GHI as fixed effect, study identity as random effect" with
  *** no covariates. That is the REGISTERED model, and it is not the model whose
  *** results are reported here or passed to the LASSO stage.

Registered (unadjusted) model, for comparison:
  Pass FDR <= 0.25 : 8
  Pass FDR <= 0.05 : 4

--- limitations ---
1. Sequencing runs are joined to the GHI on Study_ID + Treatment_Group. The join
   is at ARM level: 05_maaslin2.R sets JOIN_LEVEL <- "arm", and the key-building
   function therefore omits Segment. Several runs share one biomarker row,
   because papers report a single group mean per arm. Every run in an arm
   carries that arm's GHI. This is outcome measurement error and it caps the
   achievable correlation independently of the microbiome signal.

   *** CORRECTION. An earlier version of this file stated the join was on
   *** Study_ID + Treatment_Group + Segment. It is not. Segment enters the model
   *** as a covariate but is not part of the join key. Setting
   *** JOIN_LEVEL <- "segment" performs the organ-matched join and yields 758
   *** runs / 58 arms / 19 studies, against 1625 runs / 142 arms / 42 studies
   *** at arm level. (The LASSO stage then drops single-arm studies, leaving
   *** 141 arms across 41 studies; the two figures describe different stages.)

2. Batch correction (MMUPHin) has already been applied, and study identity
   also enters here as a random effect. This is deliberately conservative -
   both are specified in the registered protocol - but it may remove some
   genuine between-study biological variation along with technical variation.

3. FDR 0.25 is MaAsLin2's default and is permissive for a screening step.
   The 0.05 count is reported alongside so the effect of that choice is
   visible rather than buried.

4. The covariate adjustment is a deviation from the registered protocol, which
   specifies MaAsLin2 with study identity as a random effect and no additional
   fixed effects. The adjustment is defensible on its own terms - the pipeline's
   own PERMANOVA shows age explains far more community variance than the GHI
   (R2 0.11 vs 0.002) - but it materially changes the candidate pool: 8 genera
   at FDR <= 0.25 unadjusted, 28 adjusted. The 28-genus pool that feeds
   06_lasso_selection.R comes from the ADJUSTED model. Both are reported, and
   the deviation is filed as PROSPERO amendment A16.

--- provenance ---
The results in this file are identical to r_analysis/gdi/maaslin2_age_class_seg/
(verified: GHI-term coefficients match to numerical precision). The unadjusted
registered run is in r_analysis/gdi/maaslin2/. The top-level
maaslin2_all_results.tsv consumed by the LASSO stage is the adjusted run.
PDIEOF
echo "written: $D/maaslin2_summary.txt"
wc -l "$D/maaslin2_summary.txt"
