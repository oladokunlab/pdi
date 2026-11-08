# Sequence processing and analysis pipeline

**Corpus:** 50 studies, 3,037 sequencing runs, 2,064 after bioinformatic QC
**PROSPERO:** CRD420261384137
**QIIME 2:** 2024.10.1 · **SILVA:** 138 (99% OTUs, full-length, SHA256-verified)
**R:** 4.5.2

Set `PDI` to the repository root before running anything:

```bash
export PDI=/path/to/pdi
cd "$PDI/scripts"
```

---

## Preflight checks

Run these before every fresh deployment. If any check fails the working copy is
incomplete; do not proceed.

```bash
awk -F'\t' 'NR>1{n++; s+=$2} END{print n" studies, "s" runs"}' studies.tsv
#   expect: 50 studies, 3037 runs

grep "SBATCH --array" 02_denoise_per_study.sh     # expect 1-50%8
test -f 03b_build_phylogeny.sh && echo 'phylogeny script present'
grep -c "forward_only" 01_build_manifests.sh      # expect >= 2
grep -c "candidate_ft" R_scripts/00_r_config.R    # expect >= 2
```

---

## Environment

Two constraints, both of which cause silent failures if ignored.

**1. Do not set `R_LIBS_USER` in a shell that runs QIIME 2 or `sbatch`.** The R
library path leaks into the QIIME 2 environment and breaks its Python bindings.
Use separate terminals:

```bash
# QIIME 2 terminal, no R_LIBS_USER
module load QIIME2/2026.1-Amplicon

# R terminal, separate shell
module load GCC/14.3.0 R/4.5.2
export R_LIBS_USER=/path/to/R/library
```

**2. Load R by explicit version.** A bare `module load R` resolves to an older
build on some systems and the analysis scripts will fail on missing functions.

Compute nodes have no outbound network access. Install R packages from a login
node before submitting jobs.

---

## Run order

| Stage | Script | Typical runtime |
|---|---|---|
| 0 configuration | `00_config.sh` (sourced, not run) | n/a |
| 0 preflight | `00_preflight_check.sh` | seconds |
| 1 manifests | `01_build_manifests.sh` | minutes |
| 2 denoising | `02_denoise_per_study.sh` (array, 8 concurrent) | 1-3 h |
| 2 denoise summary | `02b_denoise_summary.sh` | seconds |
| 3 merge and taxonomy | `03_merge_and_taxonomy.sh` | 3-4 h |
| 4 filtering and export | `04_filter_and_export.sh` | ~10 min |
| 3a ordination | `03a_merge_and_pcoa.sh` (optional) | ~30 min |
| 3b phylogeny | `03b_build_phylogeny.sh` (optional) | 10-20 h |
| 3c refresh PERMANOVA | `03c_refresh_pcoa_permanova.sh` (optional) | ~30 min |
| 5 functional prediction | `05_picrust2.sh` (optional) | 1-2 h |
| 6 R analysis | `R_scripts/` in numerical order | minutes to hours |
| 7 verification | `run_all.sbatch` | ~20 min |

Stage 2 runs one study per array task, so a single failing study does not stop
the others. Re-run only the failed task IDs.

Retrieval helpers used while assembling the corpus are kept alongside the main
stages and are not part of a clean run: `01a_normalize_raw_data.sh`,
`01b_flatten_ngdc.sh`, `01c_reconcile_counts.sh`, `01e_reconcile_accessions.sh`,
`01f_strip_sra_suffix.sh`, `01g_subsample_wassie.sh`, and `01d_primer_audit.py`.

---

## Folder layout

| Path | Contents |
|---|---|
| `scripts/` | the bash pipeline, stages 0 to 7, plus release tooling |
| `scripts/R_scripts/` | the numbered analysis pipeline and nothing else |
| `scripts/R_scripts/data_corrections/` | corrections to extracted data, each verified against the source publication |
| `scripts/R_scripts/adhoc_checks/` | one-off audits written during development; not called by any numbered script |
| `scripts/R_scripts/superseded/` | scripts that no longer run, retained because a deposited output still cites them |
| `scripts/R_scripts/v2_and_fusion/` | the SILVA 138.2 V3-V4 sensitivity re-classification and the block-fusion analysis |
| `scripts/v2/` | the bash side of that sensitivity re-classification |

Within `R_scripts/`, the numbering groups the work: `00` to `03` build the genus
table and remove per-study batch effects; the `04` series constructs and
validates the gut health index; `05` to `48` cover the dysbiosis index,
community-level tests, and functional prediction; and the `50` series is the
verification and sensitivity suite, including risk of bias and small-study
checks.

Every R script resolves its paths from the repository root, so run them from
`R_scripts/` and they will write to `r_analysis/` at the top level.

---

## Expected warnings

These appear in normal runs and do not indicate a problem:

- `Paired-end reads not merging` for single-end studies, which are processed
  forward-only, controlled by `forward_only` in the manifest builder.
- Studies reported as `WARN` at the primer-resolution step where the publication
  does not state a primer pair. Resolution falls back to the amplicon region
  recorded in `studies.tsv`.
- Non-zero exit from the optional phylogeny stage on very large feature tables;
  the tree is not required by any reported analysis.
- `19_sparse_pls.R` will not run on a system without mixOmics. Its dependency
  RSpectra requires a compiler flag the available GCC does not provide, so
  `19b_sparse_pls_fallback.R`, which uses the CRAN `spls` package, produced the
  deposited result under `r_analysis/gdi/spls_fallback/`.

---

## If the corpus changes

Edit `studies.tsv`, then re-run from stage 3 onward. Stages 1 and 2 are
per-study, so existing studies do not need reprocessing; build manifests and
denoise only the new ones, then merge, classify, filter, and export the whole
set together. Every R analysis downstream of the genus table must then be
re-run, because batch correction is fitted across all studies at once and its
output changes when the set of studies changes.
