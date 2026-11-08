#!/bin/bash
#SBATCH --job-name=fastqc_pdi
#SBATCH --output=/scratch/user/%u/pdi/qc_reports/fastqc_%j.out
#SBATCH --error=/scratch/user/%u/pdi/qc_reports/fastqc_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=medium

# =============================================================================
# run_fastqc_multiqc.sh
#
# Quality-checks all raw FASTQ files (5,202 files, 50 studies) with FastQC,
# then rolls every report into one MultiQC summary.
#
# SUBMIT:  sbatch run_fastqc_multiqc.sh
# CHECK:   squeue -u $USER
# LOGS:    /scratch/user/$USER/pdi/qc_reports/fastqc_<jobid>.out / .err
#
# 16 cores / 4 hours should comfortably cover 5,202 files - FastQC processes
# roughly 30-60 files/minute/core on typical amplicon read lengths. If the job
# times out, resubmit with --time=08:00:00 and check the .out log to see how
# far it got; FastQC skips files that already have a report, so a resubmit is
# safe and will not redo completed work.
# =============================================================================
set -euo pipefail

module purge
module load GCC/13.2.0 OpenMPI/4.1.6 MultiQC/1.27.1
module load FastQC/0.12.1-Java-11

RAW="/scratch/user/$USER/pdi/raw_data"
OUT="/scratch/user/$USER/pdi/qc_reports"
FASTQC_OUT="${OUT}/fastqc"

mkdir -p "$FASTQC_OUT"

echo "[$(date)] Listing FASTQ files"
find "$RAW" -iname "*.fastq.gz" > "${OUT}/fastq_list.txt"
N=$(wc -l < "${OUT}/fastq_list.txt")
echo "  ${N} files found"

echo "[$(date)] Running FastQC (${SLURM_CPUS_PER_TASK} threads)"
fastqc -o "$FASTQC_OUT" -t "${SLURM_CPUS_PER_TASK}" $(cat "${OUT}/fastq_list.txt")

echo "[$(date)] Running MultiQC"
multiqc "$FASTQC_OUT" -o "$OUT" -n multiqc_report --force

echo "[$(date)] Done."
echo "Summary report: ${OUT}/multiqc_report.html"
