#!/usr/bin/env bash
# =============================================================================
# run_multiqc_per_study.sh
#
# The pooled multiqc_report.html averages quality metrics across all 5,202
# files from 50 studies at once - four different amplicon regions (V3-V4, V4,
# V3-V5, and one full-length PacBio study), different platforms, different
# labs. That pooling can obscure real differences: a PacBio read's length and
# error profile looks nothing like an Illumina short read, and mixing them
# into one quality-by-position plot is not informative for either.
#
# This builds one MultiQC report PER STUDY from the same underlying FastQC
# output already generated, so R1-vs-R2 and quality-by-position can be
# checked within a consistent platform/chemistry, not averaged across four.
#
# RUN:  bash run_multiqc_per_study.sh
# NEEDS: qc_reports/fastqc/ to already contain FastQC output (from
#        run_fastqc_multiqc.sh)
# =============================================================================
set -uo pipefail

RAW="/scratch/user/$USER/pdi/raw_data"
FASTQC_DIR="/scratch/user/$USER/pdi/qc_reports/fastqc"
OUT_DIR="/scratch/user/$USER/pdi/qc_reports/per_study"

if [ ! -d "$FASTQC_DIR" ] || [ -z "$(ls -A "$FASTQC_DIR" 2>/dev/null)" ]; then
  echo "No FastQC output found in $FASTQC_DIR"
  echo "Run run_fastqc_multiqc.sh first."
  exit 1
fi

mkdir -p "$OUT_DIR"
module load GCC/13.2.0 OpenMPI/4.1.6 MultiQC/1.27.1 2>/dev/null || true

n_studies=0
for study_dir in "$RAW"/*/; do
  study=$(basename "$study_dir")
  # FastQC output filenames are derived from the FASTQ filenames, which are
  # SRR-accession-based, not study-name-based. Find which FastQC reports
  # belong to this study by matching accession prefixes present in this
  # study's raw_data folder.
  accessions=$(find "$study_dir" -iname "*.fastq.gz" -exec basename {} \; \
               | sed -E 's/(_[12])?\.fastq\.gz$//' | sort -u)
  n_acc=$(echo "$accessions" | grep -c . || true)
  if [ "$n_acc" -eq 0 ]; then continue; fi

  study_out="${OUT_DIR}/${study}"
  mkdir -p "$study_out"

  # Symlink (not copy) the matching FastQC reports into a per-study staging
  # folder so MultiQC scans only this study's files.
  stage="${study_out}/_fastqc_staged"
  mkdir -p "$stage"
  for acc in $accessions; do
    find "$FASTQC_DIR" -iname "${acc}*fastqc*" -exec ln -sf {} "$stage/" \;
  done

  if [ -z "$(ls -A "$stage" 2>/dev/null)" ]; then
    echo "  ${study}: no matching FastQC reports found, skipping"
    rm -rf "$stage"
    continue
  fi

  multiqc "$stage" -o "$study_out" -n "${study}_multiqc" --force --quiet
  rm -rf "$stage"
  n_studies=$((n_studies + 1))
  echo "  ${study}: report written (${n_acc} accessions)"
done

echo ""
echo "=== DONE ==="
echo "${n_studies} per-study reports written under: ${OUT_DIR}/<study>/<study>_multiqc.html"
echo "Pooled report (all studies, use for a whole-corpus sanity check only): "
echo "  /scratch/user/$USER/pdi/qc_reports/multiqc_report.html"
