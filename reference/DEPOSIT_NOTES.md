# Deposit notes — read before using the extraction workbook

## The `*_Unit` columns are per-study templates, not per-measurement labels

`GDI_Metadata.xlsx`, sheet `Sample_with_Zscores`, carries a `<marker>_Unit`
column for most biomarkers. These record the unit convention of the STUDY and
are broadcast across every marker column for that study's rows. They do NOT
indicate which markers a study measured.

The evidence is direct: in all 80 rows that carry a unit for nine different
markers, every marker is given the identical unit string. For example every
Cao_2020 row reads `ng/L` for IL6_Unit, TNFa_Unit, IL1b_Unit, IFNg_Unit,
MUC2_Unit and ZO1_Unit alike — a unit that cannot apply to a cytokine
concentration, a mucin and a tight-junction protein simultaneously.

Two consequences:

* **Do not infer coverage from a unit column.** Twelve studies carry an
  `IFNg_Unit` with no IFN-γ value, and the same pattern appears for IL-6 (8
  studies), IL-1β (8), TNF-α (10), ZO-1 (10) and MUC2 (10). This is the
  template, not missing extraction.
* **Coverage is the value columns only.** Marker coverage for the whole panel
  is tabulated in `marker_coverage.csv`.

Unit heterogeneity does not affect the index: every marker is z-scored within
study before orientation, so units cancel. That holds only where a marker's
rows within a study share a measurement compartment, and in two studies they
did not. Zhang_2022 held serum IL-1beta and TNF-alpha in pg/mL alongside
mucosal fold-change values, and Ding_2025 held serum IgG and IgM in ug/mL
alongside ileal ug/mg protein. In both, the non-systemic rows were blanked in
the bare marker column and retained in the `*_mucosa_fc` and `*_mucosa`
columns, so no value was discarded. The systemic compartment is the intended
one: the domain is `Immunoglobulin_Systemic` and its third marker is
`IgA_serum`.

## Derived marker columns

Several markers appear in three parallel forms — a bare name (`IL6`), a
unit-converted variant (`IL6_pg_mL`) and a z-scored variant (`IL6_zscore`).
These are the same underlying measurement. **The bare names are the canonical
ones used in the index.** Any inventory must deduplicate them or it will
triple-count. Verified for IFN-γ: all 22 `IFNg_pg_mL` values already appear in
`IFNg`, and no row carries an `IFNg_conv_status` other than `no_data` without a
value in `IFNg`.

## Gut permeability

FITC-dextran and D-lactate were registered outcomes and remained eligible, but
no value for either was extracted from any included study. The `FITC` column
exists and is empty; there is no D-lactate column.

## Quarantined script

`superseded/06f_challenge_gdi.R` does not parse and has never run. See the
header block inside the file. Nothing in the review depends on it.
