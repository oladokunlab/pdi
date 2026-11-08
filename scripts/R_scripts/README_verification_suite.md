# PDI verification and robustness suite — R scripts

Everything in this directory reproduces, as plain `Rscript` files, analyses that were
originally run interactively or in Python. Nothing here depends on any interactive
environment: each script reads only from the archived pipeline outputs, prints its
results, and writes a `.txt` log plus a `.csv` of the numbers to
`$PDI_ROOT/r_analysis/verification/`.

Three scripts in the original pipeline (31, 40, 42) computed results but wrote no
output file. Every script here calls `open_log()` first, so that cannot recur.

## Running

On a Grace login node, **once**, to install packages (compute nodes have no internet):

```
tmux new -s pdi
module purge && module load GCC/14.3.0 R/4.5.2-gfbf-2025b
export R_LIBS_USER=/scratch/user/$USER/R/library && mkdir -p $R_LIBS_USER
Rscript -e 'install.packages(c("vegan","zCompositions","glmnet","lme4","lmerTest","readxl"),
            lib=Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")'
```

Detach with `Ctrl-b d`, reattach with `tmux attach -t pdi`.

Then either submit the whole suite:

```
sbatch run_all.sbatch
```

or run one script at a time:

```
export PDI_ROOT=/scratch/user/$USER/pdi
export PDI_WORKBOOK=$PDI_ROOT/scripts/GDI_Metadata.xlsx
export BACTOTRAITS_CSV=$PDI_ROOT/reference/BACTOTRAITS_DATASET_2026-01-28_GENUSLVL.csv
Rscript 61_ghi_scale_and_sensitivity.R
```

`PDI_NO_SINK=1` disables the log tee, for sourcing inside RStudio.

## Contents

| Script | What it does | Runtime |
|---|---|---|
| `50_setup.R` | shared paths, loaders, within-study helpers — sourced by all others | — |
| `51_batch_correction_check.R` | the registered batch-correction verification | heavy |
| `54_permanova_five_configs.R` | community vs GHI under five configurations, with age as control | heavy |
| `57_scfa_guild.R` | butyrate-producer guild vs measured caecal SCFA, with LOSO | fast |
| `58_fb_ratio.R` | Firmicutes:Bacteroidetes vs GHI and vs treatment class | ~3 min |
| `59_displacement_coupling.R` | replaces script 42: does a change in community track a change in health | fast |
| `60_ghi_reliability.R` | Cronbach's alpha by domain and the attenuation ceiling | fast |
| `61_ghi_scale_and_sensitivity.R` | validated contrasts, the 0–10 scale, Ig orientation sensitivity | fast |
| `63_rob_and_publication_bias.R` | registered RoB and platform sensitivity, funnel-asymmetry proxy | fast |
| `65_multiplicity.R` | Benjamini–Hochberg within each test family | fast |

The two heavy scripts each run `cmultRepl` on the full 2,064 × 115 table and need
roughly 16–48 GB; run them on a compute node, not a login node.

## Verified output

Every script carries its expected numbers in the header comment. The following were
confirmed to reproduce exactly on the archived data:

| Result | Value |
|---|---|
| Composite = mean of available domains | reproduced to 6.7 × 10⁻¹⁵ |
| Additive contrast | +0.8177, k = 36, p = 5.9 × 10⁻⁸ |
| Challenge contrast | +0.8244, k = 24, p = 0.0039 |
| 0–10 scale | ρ = 1.000 against the existing 0–100 column; 3 of 529 clipped |
| Ig orientation (+1 / −1 / excluded) | +0.818 / +0.561 / +0.785 — conclusions invariant |
| Domain alpha | Barrier 0.668, Inflammation 0.729, SCFA 0.889, Ig 0.773 |
| Displacement coupling, signed | −0.186, p = 0.0399; unsigned +0.040, p = 0.504 |
| Butyrate guild | ρ = +0.488, p = 0.013, q = 0.038; 10/10 LOSO folds positive |
| F:B ratio | ρ = −0.148, p = 0.11; contrasts +0.005 and +0.041 |
| Risk of bias | two zero-variance domains; 7 High across 7 studies |
| Sensitivity analyses | all four reproduce (see below) |

## Three numbers this suite corrected

Writing the scripts caught three errors in the interactive working notes. The script
values are the ones to report.

1. **Illumina-only sensitivity** is additive **+0.849 (k = 33)** and challenge
   **+0.826 (k = 23)** — not +0.775 / +0.869. The study sets matched; the earlier
   effect sizes did not.
2. **The performance domain** uses `BWG_g` (weight *gain*), which correlates with FCR
   at only **+0.05** — not +0.49. That figure belongs to `BW_g` (absolute weight), a
   different column. Only the morphology domain is a difference of strongly correlated
   markers (VH vs CD, +0.52).
3. **The butyrate guild** gives **ρ = +0.488, q = 0.038** under the deposited code,
   against +0.495 / 0.048 in the working notes. Same 33 arms and 10 studies, same
   conclusion; the script value is the reproducible one.

## Two things the scripts deliberately do not hide

- The registered **Egger test** is not computable: published reports give no per-arm
  standard errors. `63_` runs a precision proxy instead and says so.
- The registered **n ≥ 6 per group** sensitivity analysis is not executable: the
  extraction recorded sample size as free text. `63_` states this rather than omitting it.

## Not yet converted

`52_feature_restoration.R`, `53_stability_calibrated.R`, `55_balance_search_loso.R`,
`56_fusion.R` and the PICRUSt2 screen exist as analyses but not yet as scripts. They are
all bounding analyses for the null; the numbers are in the manuscript outline. Convert
them before deposit if a reviewer asks for the full search to be reproducible.
