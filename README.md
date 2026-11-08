# Poultry Gut Health and Dysbiosis Indices

Analysis code for a systematic review and meta-analysis of gut-health biomarkers
and 16S rRNA community data in poultry.

**Registration:** PROSPERO [CRD420261384137](https://www.crd.york.ac.uk/prospero/display_record.php?RecordID=1384137)

Two indices are constructed and evaluated:

- **GHI (Gut Health Index)** — a seven-domain composite built from biomarkers
  extracted from published trials. It separates treatment from control in both
  directions: beneficial additives raise it by +0.821 SD (k = 36) and challenges
  lower it by +0.868 SD (k = 25), pooled by random effects weighted on arm count
  (`04e_mixed_model_validation.R`). The unweighted study-mean series in
  `63_rob_and_publication_bias.R` gives +0.820 and +0.846 over the same studies;
  the two series are reported separately and never mixed.
- **GDI (Gut Dysbiosis Index)** — an attempt to predict the GHI from 16S rRNA
  community composition. It does not validate. Nine independent methods were
  applied and the result is reported as a null.

One positive sequence-level result is reported: a butyrate-producer guild
log-ratio tracks directly measured caecal butyrate (ρ = +0.489, q = 0.038,
direction retained in 10 of 10 leave-one-study-out folds) while failing to track
the health index on the same arms (ρ = +0.233, p = 0.267).

## Corpus

| | |
|---|---|
| Included studies | 50 |
| Sequencing runs deposited | 3,037 |
| Runs after bioinformatic QC | 2,064 |
| Runs in the analysed matrix | 1,625 (42 studies) |
| Samples scored on the GHI | 529 |

`k` counts independent experiments, not publications. Ma_2024 reports two
cohorts from one paper and one BioProject (chicks on a faecal-transplant trial,
hens on an alpha-linolenic-acid trial), carried as `Ma_2024` and
`Ma_2024_Layers`, so additive k = 36 arises from 35 publications and challenge
k = 25 from 24. The challenge contrast additionally excludes Massacci_2019,
whose control arms precede its challenge rather than running alongside them.

Two properties of the deposited matrices are load-bearing and easy to miss.
The 84-genus CLR matrix is effectively caecal: the 277 samples lost to
zero-replacement are strongly non-random by intestinal segment, so community
results should be read as caecal rather than whole-gut. And the PICRUSt2
predictions carry a weighted NSTI of 0.057 on average with a maximum of 0.185,
which exceeds the 0.15 bound PICRUSt2's own documentation gives for reliable
inference; the functional-prediction check is reported as supporting evidence
only for that reason.

No raw sequencing data is distributed here. Every included study has a public
accession; these are listed in `scripts/studies.tsv`.

## Layout

    scripts/                  shell pipeline (QIIME 2), configs, study manifests
    scripts/R_scripts/        analysis and verification scripts
    data/                     input data (see below)
    reference/                external reference data, the index definition,
                              worked examples, and notes on reading the workbook
    results/                  reference outputs from the verification suite,
                              for comparison against a fresh run

`data/` holds four files:

| File | What it is |
|---|---|
| `GDI_Metadata.xlsx` | the extraction workbook: biomarker values, study and arm metadata, risk-of-bias assessments |
| `GDI_metadata_qiime2.tsv` | per-sequencing-run sample metadata |
| `studies.tsv` | the 50 included studies with their public accessions |
| `genus_counts.tsv` | genus-level count table, the entry point for the statistical analysis |

Raw sequencing data is not distributed here. Every included study has a public
accession, listed in `data/studies.tsv`.

`genus_counts.tsv` is the output of the QIIME 2 stages and of
`01_genus_agglomeration.R`. Shipping it means the statistical analysis — batch
correction, index construction, and everything downstream — reproduces from this
repository alone, without re-running sequence processing. The shell pipeline and
`01_genus_agglomeration.R` are included for completeness and can be run from the
public accessions if you want to rebuild that table.

### What is not in this repository

The data-correction scripts are excluded. Reproducing the published results
requires the corrected workbook, which is in `data/`, not the history of how it
was corrected. `DATA_CORRECTIONS.md` records what was changed, in what order,
and the state of each correction against the distributed workbook. The scripts
themselves, and timestamped copies of the workbook before each change, are in
the full archive deposited alongside the publication.

Exploratory scripts whose results are not reported are also excluded.

## Requirements

- QIIME 2 2024.10.1, SILVA 138 (99% OTUs, full-length)
- R 4.5.2 with `data.table`, `vegan`, `zCompositions`, `glmnet`, `lme4`,
  `lmerTest`, `readxl`, `Maaslin2`, `MMUPHin`, `metafor`, `ranger`, `mixOmics`,
  `spls`. The last four are hard dependencies of result-producing scripts, not
  optional: `metafor` in 18, 26 and 33, `ranger` in 25, `mixOmics` in 19 and
  `spls` in 19b. `rstudioapi` is genuinely optional and guarded.
- PICRUSt2 (optional; used for the functional-prediction check only)

## Running the analysis

The sequence-processing stages require a cluster and roughly 20 CPU-hours. The R
stages run on a workstation in minutes.

```bash
# sequence processing, in order
bash scripts/00_preflight_check.sh
sbatch scripts/01_build_manifests.sh
sbatch scripts/02_denoise_per_study.sh      # array job, 1-3 h
sbatch scripts/03_merge_and_taxonomy.sh     # 3-4 h
bash  scripts/04_filter_and_export.sh

# R analysis
Rscript scripts/R_scripts/01_genus_agglomeration.R
Rscript scripts/R_scripts/03_batch_correction.R
Rscript scripts/R_scripts/04b_ghi_pca.R
Rscript scripts/R_scripts/05_maaslin2.R
Rscript scripts/R_scripts/06_lasso_selection.R
Rscript scripts/R_scripts/08_loso_validate.R

# independent verification of every reported number
sbatch scripts/run_all.sbatch
```

The six R lines above are the construction path: they build the genus table,
correct for study batch, fit the index and run the three selection analyses.
They are not the whole analysis. Roughly forty further scripts in
`scripts/R_scripts/` produce deposited results -- the sensitivity, diagnostic,
functional and trait-based analyses among them -- and are run together by
`run_all.sbatch`.

`scripts/README.md` gives the full stage-by-stage guide, expected runtimes and
the preflight checks. `scripts/R_scripts/README_verification_suite.md` documents
the verification suite. Two environment switches do not default to the
deposited configuration and must be set explicitly to reproduce it:
`CLDN1_MODE=merge` for `04b_ghi_pca.R` (its default drops Cldn1 and gives a
17-marker construction) and `AGE_COVARIATE=no` for `03_batch_correction.R` (its
default turns age on, which confounds with study and halts MMUPHin).

## Scoring a new study on the GHI

The index is portable. `reference/ghi_index_definition.tsv` carries the marker
panel, the direction of effect for each marker, the PC1 loadings and the
calibration constants; `scripts/R_scripts/score_new_ghi.R` reads it at run time.

```bash
Rscript scripts/R_scripts/score_new_ghi.R my_biomarkers.tsv
```

Input is a tab-separated table with a `Study_ID` column, one row per treatment
group including the control arm, and any subset of the 18 panel markers as
numeric columns. A worked example is in `reference/example_ghi_input.tsv`.

**Interpretation.** Markers are standardised within study, so scores are relative
to the baseline of the data supplied, not to an external reference population.
Higher is better; zero is the mean of the supplied data. Scores from separate
runs are not comparable, and no diagnostic cut-point is proposed — see
`reference/DEPOSIT_NOTES.md`.

## Notes before using the extraction workbook

`reference/DEPOSIT_NOTES.md` documents three things that are easy to misread: the
`*_Unit` columns are per-study templates rather than per-measurement labels,
several markers appear in parallel raw/converted/z-scored forms that must be
deduplicated, and the registered gut-permeability markers carry no extracted
values.

## Citation

Please cite the accompanying publication. The analysis outputs are deposited
separately with a dataset DOI; raw sequencing data is not redeposited, as every
included study has a public accession listed in `data/studies.tsv`.

## Licence

Code is released under the MIT Licence (`LICENSE`). Extracted biomarker values
remain the property of the original publications and are provided for
verification of the analyses reported here.
