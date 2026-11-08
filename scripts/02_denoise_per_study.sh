#!/bin/bash
#SBATCH --job-name=pdi_denoise
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-50%8
#SBATCH --output=/scratch/user/%u/pdi/logs/denoise_%A_%a.log
# =============================================================================
# 02_denoise_per_study.sh
#
# Imports raw reads, trims primers, and denoises one study per array task.
#
# Usage:
#   sbatch 02_denoise_per_study.sh
#
# Each study is processed as an independent SLURM array task rather than in
# one serial loop, so that a failure or timeout on the largest study does
# not affect the others, and studies of very different sizes do not compete
# for the same wall-clock budget. The %8 suffix on the array declaration
# limits concurrent tasks to 8.
#
# Denoising is performed per study rather than on a single pooled run.
# DADA2 learns a sequencing error model from the data it is given, and that
# model is only valid within one sequencing run/chemistry. The corpus spans
# multiple 16S regions and multiple sequencing platforms; pooling raw reads
# across these would train a single error model that fits none of them.
# Studies are instead denoised independently and combined afterward at the
# feature-table level (stage 03), which is valid because ASVs are exact
# sequences and are therefore directly comparable across studies.
#
# This script is resumable: each stage checks for its own expected output
# and is skipped if already present, so the array can be resubmitted after
# a partial failure without repeating completed work.
# =============================================================================

# Under sbatch, SLURM copies this script to a spool directory and executes
# it from there, so $(dirname "$0") does not point at the scripts directory
# and a plain `source 00_config.sh` would fail. SLURM_SUBMIT_DIR is the
# directory sbatch was invoked from, which is where the repository actually
# lives; this falls back to the script's own directory for interactive use
# on the login node.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
if [[ ! -f "${SCRIPT_DIR}/00_config.sh" ]]; then
  echo "FATAL: 00_config.sh not found in ${SCRIPT_DIR}" >&2
  echo "Submit with:  cd /scratch/user/$USER/pdi/scripts && sbatch $(basename "$0")" >&2
  exit 1
fi
source "${SCRIPT_DIR}/00_config.sh"

# If sourcing succeeded but exported nothing (e.g. a truncated or corrupted
# config file), every downstream command would run against empty variables.
# Failing here avoids wasting compute allocation on a job that cannot do
# anything useful.
if [[ -z "${PDI_ROOT:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh loaded but PDI_ROOT/log() are unset." >&2
  exit 1
fi
activate_qiime

# ---- Select this task's study from studies.tsv ------------------------------
# +1 because line 1 of studies.tsv is the header. The table is sorted
# smallest-study-first, so low array indices complete quickly and surface
# any systematic failure before the largest studies have run.
LINE=$((SLURM_ARRAY_TASK_ID + 1))
# `|| true` guards against `set -e`: `read` returns non-zero when the final
# line of studies.tsv has no trailing newline, which would otherwise abort
# the task before STUDY is even checked.
IFS=$'\t' read -r STUDY N_RUNS REGION PLATFORM LAYOUT PRIMER_F PRIMER_R STATUS \
  PRIMERS_PRESENT READ_ORIENTATION \
  < <(sed -n "${LINE}p" "$STUDY_TABLE" | tr -d '\r') || true
PRIMERS_PRESENT="${PRIMERS_PRESENT:-YES}"
READ_ORIENTATION="${READ_ORIENTATION:-normal}"

[[ -z "${STUDY:-}" ]] && { log "No study on line ${LINE}; nothing to do."; exit 0; }

log "=== Task ${SLURM_ARRAY_TASK_ID}: ${STUDY} (${N_RUNS} runs, ${REGION}, ${LAYOUT}) ==="

# ---- Skip studies flagged on hold -------------------------------------------
# A study is held when a required input (typically a verified primer
# sequence) is not yet available. Processing is skipped explicitly rather
# than proceeding with an unverified assumption.
if [[ "$STATUS" == HOLD_* ]]; then
  log "HOLD: ${STUDY} is parked (${STATUS}). Not processing. See README section 'Studies on hold'."
  exit 0
fi

SDIR="${WORK_DIR}/${STUDY}"
mkdir -p "$SDIR"
MAN="${MANIFEST_DIR}/${STUDY}_manifest.tsv"
[[ -s "$MAN" ]] || { log "ERROR: manifest missing for ${STUDY}. Run 01 first."; exit 1; }

# Resumability gate
already_done "${SDIR}/table.qza" && { log "${STUDY} already denoised."; exit 0; }

# =============================================================================
# MIXED-PRIMER STUDIES (optional)
#
# Some studies contain a mix of samples where primers are still present and
# samples where the submitting laboratory already removed them. A single
# study-level primers_present flag cannot express that, and running cutadapt
# with --p-discard-untrimmed across the whole study would discard every read
# from the already-trimmed subset (no primer left to find).
#
# To handle this, a study may have a companion file:
#
#   manifests/<Study_ID>_noprimer_samples.txt
#
# containing one sample-id per line for the samples whose primers are
# ALREADY REMOVED. When that file exists, this block:
#   1. splits the study's manifest into two sub-manifests
#   2. imports and processes each with the correct primer handling
#   3. denoises each sub-batch separately
#   4. merges the resulting tables and rep-seqs back into one study output
#
# Denoising the two subsets separately is not a compromise - DADA2 learns an
# error model per run, and reads that have been through different upstream
# processing legitimately have different error profiles. Merging afterward at
# the feature-table level is valid for exactly the same reason merging across
# studies is valid: ASV identity is the exact sequence.
#
# The study keeps ONE Study_ID and ONE output directory throughout, so
# nothing downstream (metadata, batch correction, grouping) sees a split.
# =============================================================================
NOPRIMER_LIST="${MANIFEST_DIR}/${STUDY}_noprimer_samples.txt"

if [[ -s "$NOPRIMER_LIST" && "$LAYOUT" == "paired" ]]; then
  log "MIXED-PRIMER MODE: $(wc -l < "$NOPRIMER_LIST") sample(s) listed as already primer-trimmed"

  if already_done "${SDIR}/table.qza"; then
    log "${STUDY} already denoised."
    exit 0
  fi

  MAN_WITH="${SDIR}/manifest_withprimer.tsv"
  MAN_WITHOUT="${SDIR}/manifest_noprimer.tsv"

  head -1 "$MAN" > "$MAN_WITH"
  head -1 "$MAN" > "$MAN_WITHOUT"
  # grep -F -x -f matches whole lines literally, so a sample-id cannot
  # partially match a longer one (e.g. RA1 vs RA10).
  awk -F'\t' 'NR>1' "$MAN" | while IFS=$'\t' read -r sid f r; do
    if grep -qFx "$sid" "$NOPRIMER_LIST"; then
      printf '%s\t%s\t%s\n' "$sid" "$f" "$r" >> "$MAN_WITHOUT"
    else
      printf '%s\t%s\t%s\n' "$sid" "$f" "$r" >> "$MAN_WITH"
    fi
  done

  N_WITH=$(( $(wc -l < "$MAN_WITH") - 1 ))
  N_WITHOUT=$(( $(wc -l < "$MAN_WITHOUT") - 1 ))
  log "  with primers: ${N_WITH} sample(s)   already trimmed: ${N_WITHOUT} sample(s)"
  if [[ "$N_WITH" -eq 0 || "$N_WITHOUT" -eq 0 ]]; then
    log "ERROR: mixed-primer mode requires both subsets to be non-empty."
    log "       If ALL samples are pre-trimmed, set primers_present=NO in studies.tsv instead."
    exit 1
  fi

  SUB_TABLES=(); SUB_SEQS=()
  for subset in with without; do
    if [[ "$subset" == "with" ]]; then
      SUBMAN="$MAN_WITH"; TAG="withprimer"
    else
      SUBMAN="$MAN_WITHOUT"; TAG="noprimer"
    fi

    $QIIME tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path "$SUBMAN" \
      --output-path "${SDIR}/reads_${TAG}.qza" \
      --input-format PairedEndFastqManifestPhred33V2

    if [[ "$subset" == "with" ]]; then
      log "  [${TAG}] cutadapt trim-paired  F=${PRIMER_F}  R=${PRIMER_R}"
      $QIIME cutadapt trim-paired \
        --i-demultiplexed-sequences "${SDIR}/reads_${TAG}.qza" \
        --p-cores "$NCORES" \
        --p-front-f "$PRIMER_F" \
        --p-front-r "$PRIMER_R" \
        --p-discard-untrimmed \
        --p-no-indels \
        --o-trimmed-sequences "${SDIR}/reads_${TAG}_trimmed.qza"
    else
      log "  [${TAG}] skipping cutadapt - primers already removed by submitter"
      cp "${SDIR}/reads_${TAG}.qza" "${SDIR}/reads_${TAG}_trimmed.qza"
    fi

    log "  [${TAG}] DADA2 denoise-paired"
    $QIIME dada2 denoise-paired \
      --i-demultiplexed-seqs "${SDIR}/reads_${TAG}_trimmed.qza" \
      --p-trim-left-f "$TRIM_LEFT" --p-trim-left-r "$TRIM_LEFT" \
      --p-trunc-len-f "$TRUNC_LEN" --p-trunc-len-r "$TRUNC_LEN" \
      --p-max-ee-f "$MAX_EE" --p-max-ee-r "$MAX_EE" \
      --p-chimera-method "$CHIMERA_METHOD" \
      --p-n-threads "$NCORES" \
      --o-table "${SDIR}/table_${TAG}.qza" \
      --o-representative-sequences "${SDIR}/rep-seqs_${TAG}.qza" \
      --o-denoising-stats "${SDIR}/denoising-stats_${TAG}.qza" \
      --verbose

    SUB_TABLES+=(--i-tables "${SDIR}/table_${TAG}.qza")
    SUB_SEQS+=(--i-data "${SDIR}/rep-seqs_${TAG}.qza")
  done

  log "  merging the two sub-batches into one study-level output"
  $QIIME feature-table merge "${SUB_TABLES[@]}" \
    --p-overlap-method error_on_overlapping_sample \
    --o-merged-table "${SDIR}/table.qza"
  $QIIME feature-table merge-seqs "${SUB_SEQS[@]}" \
    --o-merged-data "${SDIR}/rep-seqs.qza"

  # Denoising stats are kept per sub-batch rather than merged: the two
  # subsets went through different trimming, so a single combined stats
  # table would misrepresent the filtering step for one of them.
  log "  denoising stats kept separately: denoising-stats_withprimer.qza, denoising-stats_noprimer.qza"
  $QIIME metadata tabulate \
    --m-input-file "${SDIR}/denoising-stats_withprimer.qza" \
    --o-visualization "${SDIR}/denoising-stats_withprimer.qzv"
  $QIIME metadata tabulate \
    --m-input-file "${SDIR}/denoising-stats_noprimer.qza" \
    --o-visualization "${SDIR}/denoising-stats_noprimer.qzv"

  log "=== ${STUDY} COMPLETE (mixed-primer mode) ==="
  exit 0
fi

# =============================================================================
# STEP 1 - IMPORT
# Converts the manifest-referenced FASTQ files into a QIIME 2 artifact
# (.qza). A .qza retains provenance metadata recording every command that
# produced it, supporting later audit of the processing history.
# =============================================================================
if ! already_done "${SDIR}/reads.qza"; then
  log "[1/5] Importing ${STUDY}"
  if [[ "$LAYOUT" == "paired" ]]; then
    $QIIME tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path "$MAN" \
      --output-path "${SDIR}/reads.qza" \
      --input-format PairedEndFastqManifestPhred33V2
  else
    # ccs, single_merged, and forward_only all import as single-end - one
    # file per sample.
    $QIIME tools import \
      --type 'SampleData[SequencesWithQuality]' \
      --input-path "$MAN" \
      --output-path "${SDIR}/reads.qza" \
      --input-format SingleEndFastqManifestPhred33V2
  fi
fi

# =============================================================================
# STEP 2 - SUMMARISE RAW READS
# Produces an interactive quality-score visualization (view at
# view.qiime2.org). This is a sanity check on the import step - confirming
# expected read length and depth - and is not used to select a truncation
# length, since the registered setting (--p-trunc-len 0) applies no fixed
# truncation.
# =============================================================================
if ! already_done "${SDIR}/reads_summary.qzv"; then
  log "[2/5] Summarising raw reads"
  $QIIME demux summarize \
    --i-data "${SDIR}/reads.qza" \
    --o-visualization "${SDIR}/reads_summary.qzv"
fi

# =============================================================================
# STEP 3 - TRIM PRIMERS WITH CUTADAPT
# Primer sequences are synthetic and not part of the biological signal.
# Leaving them untrimmed distorts DADA2's error model (invariant primer
# bases are counted as informative sequence) and inflates the ASV count
# (degenerate primer positions are read as genuine sequence variation).
#
# Primer sequences are read per study from studies.tsv, since the corpus
# spans multiple 16S regions rather than a single primer pair.
#
# --p-discard-untrimmed removes any read in which the primer sequence is
# not found, which also serves as a validation check: a study losing most
# of its reads at this step indicates an incorrect primer assignment.
# =============================================================================
if [[ "$PRIMERS_PRESENT" == "NO" ]]; then
  # Some studies were deposited with primers already removed by the
  # submitting laboratory. For these, cutadapt is skipped entirely: running
  # --p-discard-untrimmed against reads with no primer sequence present
  # would discard every read and yield an empty feature table.
  log "[3/5] SKIPPING cutadapt - primers already removed by submitter"
  if ! already_done "${SDIR}/reads_trimmed.qza"; then
    cp "${SDIR}/reads.qza" "${SDIR}/reads_trimmed.qza"
    $QIIME demux summarize \
      --i-data "${SDIR}/reads_trimmed.qza" \
      --o-visualization "${SDIR}/reads_trimmed_summary.qzv"
  fi

elif [[ "$LAYOUT" == "paired" ]]; then
  if ! already_done "${SDIR}/reads_trimmed.qza"; then
    # PRIMERS_PRESENT=TAIL: for studies where the submitting laboratory
    # trimmed partway into the primer sequence, only the surviving tail is
    # present at the read start. studies.tsv holds that shortened sequence
    # as primer_f for these studies, so cutadapt is given the correct
    # (shorter) target directly.
    log "[3/5] cutadapt trim-paired  F=${PRIMER_F}  R=${PRIMER_R}  (primers=${PRIMERS_PRESENT})"
    $QIIME cutadapt trim-paired \
      --i-demultiplexed-sequences "${SDIR}/reads.qza" \
      --p-cores "$NCORES" \
      --p-front-f "$PRIMER_F" \
      --p-front-r "$PRIMER_R" \
      --p-discard-untrimmed \
      --p-no-indels \
      --o-trimmed-sequences "${SDIR}/reads_trimmed.qza"

    $QIIME demux summarize \
      --i-data "${SDIR}/reads_trimmed.qza" \
      --o-visualization "${SDIR}/reads_trimmed_summary.qzv"
  fi

elif [[ "$LAYOUT" == "forward_only" ]]; then
  # Studies where forward and reverse reads do not overlap sufficiently to
  # be joined. Only the forward read is processed, so only the forward
  # primer is trimmed; no reverse-complement adapter is needed, since the
  # forward read does not extend far enough into the amplicon to reach the
  # reverse primer binding site.
  if ! already_done "${SDIR}/reads_trimmed.qza"; then
    log "[3/5] cutadapt trim-single (forward read only)  front=${PRIMER_F}"
    $QIIME cutadapt trim-single \
      --i-demultiplexed-sequences "${SDIR}/reads.qza" \
      --p-cores "$NCORES" \
      --p-front "$PRIMER_F" \
      --p-discard-untrimmed \
      --o-trimmed-sequences "${SDIR}/reads_trimmed.qza"

    $QIIME demux summarize \
      --i-data "${SDIR}/reads_trimmed.qza" \
      --o-visualization "${SDIR}/reads_trimmed_summary.qzv"
  fi

elif [[ "$LAYOUT" == "single_merged" ]]; then
  # Studies deposited with forward and reverse reads already joined into a
  # single read per sample (approximately 400-470 bp, spanning the full
  # amplicon). A merged read has the structure [forward primer][insert]
  # [reverse complement of reverse primer], so both ends require trimming:
  # --p-front removes the forward primer, and --p-adapter removes the
  # reverse complement of the reverse primer from the 3' end. Leaving
  # either untrimmed would distort DADA2's error model in the same way an
  # untrimmed primer would on paired-end data.
  if ! already_done "${SDIR}/reads_trimmed.qza"; then
    REVCOMP_R=$(python3 -c "
comp = str.maketrans('ACGTMRWSYKVHDBNacgtmrwsykvhdbn','TGCAKYWSRMBDHVNtgcakywsrmbdhvn')
print('${PRIMER_R}'.translate(comp)[::-1])
")
    log "[3/5] cutadapt trim-single (merged read)  front=${PRIMER_F}  adapter(revcomp R)=${REVCOMP_R}"
    $QIIME cutadapt trim-single \
      --i-demultiplexed-sequences "${SDIR}/reads.qza" \
      --p-cores "$NCORES" \
      --p-front "$PRIMER_F" \
      --p-adapter "$REVCOMP_R" \
      --p-discard-untrimmed \
      --o-trimmed-sequences "${SDIR}/reads_trimmed.qza"

    $QIIME demux summarize \
      --i-data "${SDIR}/reads_trimmed.qza" \
      --o-visualization "${SDIR}/reads_trimmed_summary.qzv"
  fi
fi
# PacBio CCS studies (layout=ccs) do not run cutadapt as a separate step;
# denoise-ccs accepts primer sequences directly via --p-front/--p-adapter
# and performs read orientation internally.

# =============================================================================
# STEP 4 - DENOISE WITH DADA2
#
# DADA2 filters reads by expected error, learns a per-study error model,
# dereplicates reads, infers exact amplicon sequence variants (ASVs), joins
# forward and reverse reads where applicable, and removes chimeric
# sequences.
#
# ASVs are used in preference to OTU clustering because ASV identity is
# determined by the exact sequence rather than by clustering relative to a
# particular dataset. This is what permits studies to be denoised
# independently and then combined at the feature-table level in stage 03:
# an identical ASV sequence produced in two different studies is
# recognized as the same feature.
#
# --p-trunc-len-f/r 0 applies no fixed-length truncation (the registered
# setting; see 00_config.sh). --p-max-ee-f/r sets the maximum tolerated
# expected error per read, adapting to each read's own quality profile
# rather than applying one fixed cutoff position across all reads.
#
# denoise-paired requires reads that have not been pre-merged or
# pre-quality-filtered - it performs both steps internally using the raw
# quality scores to fit its error model. Any workflow that merges or
# quality-filters reads before this step is incompatible with DADA2.
# =============================================================================
log "[4/5] DADA2 denoising (${LAYOUT})"

if [[ "$LAYOUT" == "paired" ]]; then
  $QIIME dada2 denoise-paired \
    --i-demultiplexed-seqs "${SDIR}/reads_trimmed.qza" \
    --p-trim-left-f "$TRIM_LEFT" \
    --p-trim-left-r "$TRIM_LEFT" \
    --p-trunc-len-f "$TRUNC_LEN" \
    --p-trunc-len-r "$TRUNC_LEN" \
    --p-max-ee-f "$MAX_EE" \
    --p-max-ee-r "$MAX_EE" \
    --p-chimera-method "$CHIMERA_METHOD" \
    --p-n-threads "$NCORES" \
    --o-table "${SDIR}/table.qza" \
    --o-representative-sequences "${SDIR}/rep-seqs.qza" \
    --o-denoising-stats "${SDIR}/denoising-stats.qza" \
    --verbose

elif [[ "$LAYOUT" == "ccs" ]]; then
  # PacBio HiFi full-length 16S reads. denoise-ccs applies an error model
  # tuned for CCS read accuracy and performs primer removal internally.
  # --p-min-len/--p-max-len bracket the expected length of a full-length
  # 16S amplicon.
  $QIIME dada2 denoise-ccs \
    --i-demultiplexed-seqs "${SDIR}/reads.qza" \
    --p-front "$PRIMER_F" \
    --p-adapter "$PRIMER_R" \
    --p-min-len 1000 \
    --p-max-len 1600 \
    --p-trunc-len "$TRUNC_LEN" \
    --p-max-ee "$MAX_EE" \
    --p-chimera-method "$CHIMERA_METHOD" \
    --p-n-threads "$NCORES" \
    --o-table "${SDIR}/table.qza" \
    --o-representative-sequences "${SDIR}/rep-seqs.qza" \
    --o-denoising-stats "${SDIR}/denoising-stats.qza" \
    --verbose

elif [[ "$LAYOUT" == "forward_only" ]]; then
  # denoise-single on the forward read only. No merge step is performed,
  # since the reverse read is excluded from processing entirely for these
  # studies. Registered trunc-len and max-ee settings apply as elsewhere.
  $QIIME dada2 denoise-single \
    --i-demultiplexed-seqs "${SDIR}/reads_trimmed.qza" \
    --p-trim-left "$TRIM_LEFT" \
    --p-trunc-len "$TRUNC_LEN" \
    --p-max-ee "$MAX_EE" \
    --p-chimera-method "$CHIMERA_METHOD" \
    --p-n-threads "$NCORES" \
    --o-table "${SDIR}/table.qza" \
    --o-representative-sequences "${SDIR}/rep-seqs.qza" \
    --o-denoising-stats "${SDIR}/denoising-stats.qza" \
    --verbose

elif [[ "$LAYOUT" == "single_merged" ]]; then
  # denoise-single on the pre-merged read. No further merge step is
  # required, since forward and reverse reads were already joined prior to
  # deposit. Registered trunc-len and max-ee settings apply as elsewhere.
  $QIIME dada2 denoise-single \
    --i-demultiplexed-seqs "${SDIR}/reads_trimmed.qza" \
    --p-trim-left "$TRIM_LEFT" \
    --p-trunc-len "$TRUNC_LEN" \
    --p-max-ee "$MAX_EE" \
    --p-chimera-method "$CHIMERA_METHOD" \
    --p-n-threads "$NCORES" \
    --o-table "${SDIR}/table.qza" \
    --o-representative-sequences "${SDIR}/rep-seqs.qza" \
    --o-denoising-stats "${SDIR}/denoising-stats.qza" \
    --verbose
fi

# =============================================================================
# STEP 5 - DENOISING STATS
# Tabulates, per sample: input -> filtered -> denoised -> merged (where
# applicable) -> non-chimeric read counts. This table is inspected for
# every study before proceeding to the merge stage. Substantial read loss
# at "filtered" indicates lower-than-expected quality; loss at "merged"
# indicates insufficient read overlap for the targeted amplicon length;
# loss at "non-chimeric" indicates PCR chimera formation in the original
# library.
# =============================================================================
log "[5/5] Tabulating denoising stats"
$QIIME metadata tabulate \
  --m-input-file "${SDIR}/denoising-stats.qza" \
  --o-visualization "${SDIR}/denoising-stats.qzv"

log "=== ${STUDY} COMPLETE ==="
