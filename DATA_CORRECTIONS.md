# Data corrections applied to the extraction workbook

This file records the corrections made to `GDI_Metadata.xlsx` between first
transcription from the source publications and the workbook the analysis was run
on, and the order in which they were applied.

**The scripts listed here must not be re-run.** The distributed workbook already
contains every correction that was adopted. Several of the scripts are stale
relative to it: their target rows no longer exist, or they would re-apply a
change that was later reversed. They are retained as a record of what was
changed and why, not as a pipeline. The table below gives the state of each
against the distributed workbook, established by a dry run on 23 September 2026.
Step 3 has no dry-run mode and was verified by inspecting the workbook directly;
the evidence is given in its row.

Only the full archive carries these scripts. The code release excludes them,
because reproducing the published results requires the corrected workbook, not
the correction history.

## Order and current state

| # | Script | Effect on the distributed workbook |
|---|---|---|
| 1 | `fix_ma2024_split.R` | **Stale.** Would rename 155 sequencing runs from `Ma_2024` to `Ma_2024_FMT` in the QIIME 2 metadata. The adopted state keeps them as `Ma_2024`, matching the biomarker sheet. |
| 2 | `rename_ma2024_back.R` | No-op. Targets `Ma_2024_FMT` rows in the biomarker sheet; none exist. It cannot reverse step 1, which writes to a different file. |
| 3 | `fix_label_mismatches.R` | Applied. Reconciled treatment-group labels between `Crosswalk_16S` and `Sample_with_Zscores`. This script has no dry-run mode, so it was verified against the workbook instead: all 13 of its renamed labels are present in `Crosswalk_16S` and none of the 13 original labels remains. |
| 4 | `fix_ages_final.R` | **Partly stale.** All age corrections except one are present. It would additionally set `Zhang_2026a` to 180 days, which was applied and then reversed; see step 6. |
| 5 | `fix_missing_ages.R` | Applied. `Melaku_2024` = 42 days; `Xu_2023` mapped to days 7 and 21. |
| 6 | `fix_zhang2026a_revert.R` | **Cannot run.** Reverses the `Zhang_2026a` age change of step 4, but the 22 affected rows carry no age in their `Sample_ID`, so it refuses rather than guess. The adopted values (21 and 42 days) are in the distributed workbook. |
| 7 | `fix_classes.R` | **Not adopted.** Would reclassify `Ma_2024 Ab` and `Tang_2026 AB` as beneficial additives. Both are correct as recorded; see below. |
| 8 | `fix_classes2.R` | **Stale.** Its mapping names arms that do not exist under those names (`Xu_2023 TD`/`TFRD`; `Ma_2024 LAS`/`MAS`/`HAS`). It changes nothing. |
| 9 | `fix_iga_compartment.R` | Applied. All 64 serum IgA values moved from the mucosal immunoglobulin domain to an `IgA_serum` column. |
| 10 | `fix_siga_mislabel.R` | Applied. 15 serum IgA values removed from the `sIgA` column; 55 genuinely mucosal rows across 8 studies retained. |

Two further scripts in the same directory are read-only reports and change
nothing: `label_mismatch_check.R` and `fix_arm_names.R` (which lists arm-name
discrepancies without correcting them).

## Why `fix_classes.R` was not adopted

It matches on antibiotic tokens in the arm name and does not read the study
design. Both arms it would move are correctly classified as recorded:

**`Ma_2024 Ab`** — the study's arms are `C, Ab, C_FMT, A_FMT, S, S_Ab, S_C_FMT,
S_A_FMT`, where the `S_` prefix denotes *Salmonella* challenge. `Ab` is
therefore antibiotic given without challenge, alongside two faecal-transplant
arms, and its recorded subtype is `Unspecified_Mixed_Antibiotic`. This is
microbiota depletion before transplant, not growth promotion.
`Antibiotic_Challenge` is correct.

**`Tang_2026 AB`** — the arms are `NC` (control), `EC` (challenge), then `AB`,
`MGDL` and `MGDH`, all administered under challenge. `AB` is the antibiotic
positive control in challenged birds, and `Challenged_And_Treated` is correct.

Adopting the script would have moved five rows into the additive contrast and
changed the additive estimate from +0.818 to +0.784 SD, incorrectly. Both
figures are as at that check, before the Massacci_2019 and Xu_2023
reclassifications; the additive estimate now stands at +0.821 SD.

Antibiotic arms across the whole corpus were checked for consistency: thirteen
studies contribute one, and all thirteen are classified in accordance with their
design. Nine growth-promoter arms are `Non_Disruptive_Additive`; the two
combination arms (`Wu_2021a ANT+LPS`, `Zhang_2022 ANT+LPS`) are
`Challenged_And_Treated`.

## Recovering an earlier state

Pre-correction copies of the workbook are not distributed. The workbook shipped
here already contains every correction that was adopted, and each correction is
recorded above with the rows it changed and the source that justified it, so the
change history is readable without a second copy of the file. Re-running the
correction scripts is not a supported route to an earlier state either: several
are stale against the distributed workbook, as the table above records.


## 2026-10-24  Ding_2025 Abx reclassified; Ma_2024 S_Ab subclass refined

Ding_2025 'Abx' (vancomycin, metronidazole, kanamycin, ampicillin at 0.5:1:1:1,
100 mg/kg of feed, days 1 to 21) is described by the authors as a growth
promoter run as a positive comparator beside rosmarinic acid, not as a
depletion insult. It is reclassified from Antibiotic challenge to
Non-disruptive additive, subclass Antibiotic growth promoter. Ma_2024 'S_Ab'
keeps Infectious challenge, because the antibiotic preceded the Salmonella
rather than treating it; its subclass now records the depletion.

Applied to GDI_Metadata.xlsx (sheet Sample_with_Zscores, 4 rows),
intervention_classification_FINAL.csv, GDI_metadata_qiime2.tsv and
GDI_metadata_qiime2_class.tsv (6 runs each). Scripts 04b, 04e, 29, 61, 63 and
55 were re-run; the control anchor is unchanged at -0.222050 / 0.902578 on
n = 150 and the scored counts are 236 additive, 150 control, 68 challenged,
75 challenged-and-treated.

03_batch_correction.R protects Intervention_Class as an MMUPHin covariate.
Its published output (27 September 2026, 2,064 samples, protected covariate
Intervention_Class, no age covariate) predates the eight-class vocabulary,
which entered GDI_metadata_qiime2.tsv on 22 October 2026, so this relabel has
no bearing on the batch-corrected genus table or on any GDI result derived
from it. The script as shipped cannot be re-run against the current metadata:
the Intervention_Class column now carries the eight-class scheme, in which
Antibiotic challenge is contributed by a single study (Ma_2024) and MMUPHin
therefore reports the covariate confounded with batch. A re-run requires
either the September metadata or a covariate collapsed to the four semantic
families.

### Outputs predating the eight-class relabel

An exhaustive scan of all 206 deposited non-image outputs found 18 containing a
retired label. None carries a reported result. Three (ghi_scores_pca.tsv,
verification/ghi_scores_0_10.csv, verification/55_class_ordination.txt) retain
Challenge_Type or name the retired QIIME artifacts deliberately, as documented
in 50_setup.R. Six GDI score tables (gdi/loso_predictions.tsv,
gdi_three_ways_scores.tsv, three_way_axis_scores.tsv,
displacement_coupling_data.tsv, unsupervised_scores.tsv, axis_loso_scores.tsv)
carry the class only as an annotation column that no figure or statistic reads.
Two MaAsLin2 sensitivity runs (gdi/maaslin2_class_seg, gdi/maaslin2_age_class_seg)
used Challenge_Type as a model term and are not cited; the reported MaAsLin2
result is gdi/maaslin2_all_results.tsv, which uses the continuous index. Three
confounding diagnostics dated 13 August 2026 (confounding/) tested PCoA_Group and
the superseded six-level Intervention_Class, both retired.

### Outputs predating the eight-class relabel

An exhaustive scan of all 206 deposited non-image outputs found 18 containing a
retired label. None carries a reported result. Three (ghi_scores_pca.tsv,
verification/ghi_scores_0_10.csv, verification/55_class_ordination.txt) retain
Challenge_Type or name the retired QIIME artifacts deliberately, as documented
in 50_setup.R. Six GDI score tables (gdi/loso_predictions.tsv,
gdi_three_ways_scores.tsv, three_way_axis_scores.tsv,
displacement_coupling_data.tsv, unsupervised_scores.tsv, axis_loso_scores.tsv)
carry the class only as an annotation column that no figure or statistic reads.
Two MaAsLin2 sensitivity runs (gdi/maaslin2_class_seg, gdi/maaslin2_age_class_seg)
used Challenge_Type as a model term and are not cited; the reported MaAsLin2
result is gdi/maaslin2_all_results.tsv, which uses the continuous index. Three
confounding diagnostics dated 13 August 2026 (confounding/) tested PCoA_Group and
the superseded six-level Intervention_Class, both retired.

## 2026-09-27  Final-audit corrections to GDI_Metadata.xlsx

Five independent corrections, each verified against the source publication
before being written. A timestamped copy of the workbook was taken first.

**1. Xu_2023, arm `TFRD (Fig 5)`, reclassified to challenged and treated.**
The arm's own Notes read "TD + Total flavonoids from Rhizoma Drynariae", a
thiram-induced challenge given alongside a treatment. Recorded as a dietary
challenge it would have counted as a challenged bird in the GHI contrast. The
arm is biomarker-only, since Crosswalk_16S holds just CON and TD for this
study, so no sequencing metadata changed.

**2. Five DOIs in Inclusion_Exclusion pointed to the wrong paper.** All five
were confirmed from the source PDFs. Every journal and year in the workbook
was already correct and only the DOI column was wrong; the Chen_2026a and
Chen_2026b entries had been transposed.

**3. Wu_2021b SCFA units flagged in Notes, values left as printed.** Table 4
of the source prints "(mmol/g)" with acetate at 2555.67 and total SCFA at
5905.03, which would be 5.9 mol/g and is not physically possible. Read as
ug/g the same figures give 42.6 umol/g acetate and about 85 umol/g total,
which are ordinary broiler cecal concentrations. Nothing was overridden: the
printed values and the declared unit stand, and a note records the
discrepancy.

**4. Two studies mixed measurement compartments inside one marker column.**
`04b_ghi_pca.R` reads the bare marker column, so it did not see the workbook's
separate `*_serum` and `*_mucosa` fields. Zhang_2022 pooled four serum IL-1b
and TNF-a rows in pg/mL with eight jejunal and ileal rows expressed as fold
change, and Ding_2025 pooled three serum IgG and IgM rows in ug/mL with three
ileal rows in ug/mg protein. Within-study z-scoring cannot absorb a unit
change of that size: recomputed, the mixing moved Zhang_2022's IL-1b z-scores
from a range of -1.40 to 0.87 to a range of -0.66 to 1.71. The bare column was
blanked for the non-systemic compartment only. No value was lost, since those
rows remain in `IL1b_mucosa_fc`, `TNFa_mucosa_fc`, `IgG_mucosa`, and
`IgM_mucosa`. The systemic compartment was the one retained, because the
domain is Immunoglobulin_Systemic and its third marker is IgA_serum, and
because the corpus scores inflammation on serum and plasma concentrations
rather than on mucosal expression.

**5. Final_status used six spellings for four states.** Any count done by text
matching silently missed rows. The five variant spellings, four of `Excluded`
and one of `EXCLUDED`, were all data or metadata exclusions, each confirmed
from its own Original_reason, and were normalized to the canonical spelling.

The scripts that applied these changes are retained under
`_history/code_fixes/` and must not be re-run. The distributed workbook
already contains every correction listed here.
