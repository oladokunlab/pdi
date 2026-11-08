#!/bin/bash
#SBATCH --job-name=install_mixomics
#SBATCH --output=/scratch/user/%u/pdi/qc_reports/mixomics_install_%j.out
#SBATCH --error=/scratch/user/%u/pdi/qc_reports/mixomics_install_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=short

# =============================================================================
# install_mixomics.sh
#
# Installs BH, BiocParallel and mixOmics as a detached batch job.
#
# WHY THIS EXISTS: the last two attempts died silently mid-install with no
# error message - the process was very likely killed by a dropped terminal
# connection (SSH timeout, closed window), not a real compiler or package
# problem. Running this as an sbatch job means it keeps running on a compute
# node regardless of what happens to the terminal that submitted it.
#
# SUBMIT:  sbatch install_mixomics.sh
# CHECK:   squeue -u $USER
# LOG:     /scratch/user/$USER/pdi/qc_reports/mixomics_install_<jobid>.out
# =============================================================================
set -uo pipefail

echo "[$(date)] Clearing any stale lock from a previous interrupted install"
rm -rf /scratch/user/$USER/R/library/00LOCK*

module purge
module load GCC/14.3.0

export R_LIBS_USER=/scratch/user/$USER/R/library
R=/sw/eb/sw/R/4.5.2-gfbf-2025b/bin/Rscript

echo "[$(date)] Installing BH"
$R -e 'install.packages("BH", lib=Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")'

echo "[$(date)] Checking BH installed"
$R -e 'if (!requireNamespace("BH", quietly=TRUE)) stop("BH did not install") else cat("BH OK\n")'

echo "[$(date)] Installing BiocManager if needed, then mixOmics"
$R -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", lib=Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
BiocManager::install("mixOmics", lib=Sys.getenv("R_LIBS_USER"), update=FALSE, ask=FALSE)'

echo "[$(date)] Verifying mixOmics loads"
$R -e 'library(mixOmics); cat("mixOmics OK, version:", as.character(packageVersion("mixOmics")), "\n")'

echo "[$(date)] Done."
