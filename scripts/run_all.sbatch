#!/bin/bash
#SBATCH --job-name=pdi_verify
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
# -----------------------------------------------------------------------------
# PDI verification and robustness suite.
#
#   sbatch run_all.sbatch
#
# Compute nodes on Grace have no internet, so every package must already be
# installed in $R_LIBS_USER. Install them ONCE on a login node, inside tmux:
#
#   tmux new -s pdi
#   module purge && module load GCC/14.3.0 R/4.5.2
#   export R_LIBS_USER=/scratch/user/$USER/R/library && mkdir -p $R_LIBS_USER
#   Rscript -e 'install.packages(c("vegan","zCompositions","glmnet","lme4",
#     "lmerTest","readxl"), lib=Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")'
#
# Detach with Ctrl-b d, reattach with: tmux attach -t pdi
# -----------------------------------------------------------------------------
set -euo pipefail

module purge
module load GCC/14.3.0 R/4.5.2

export R_LIBS_USER=/scratch/user/$USER/R/library
export PDI_ROOT=/scratch/user/$USER/pdi
export PDI_WORKBOOK=$PDI_ROOT/scripts/GDI_Metadata.xlsx
export BACTOTRAITS_CSV=$PDI_ROOT/reference/BACTOTRAITS_DATASET_2026-01-28_GENUSLVL.csv

cd "$PDI_ROOT/scripts/R_scripts"

for s in 51_batch_correction_check.R \
         54_permanova_five_configs.R \
         57_scfa_guild.R \
         58_fb_ratio.R \
         59_displacement_coupling.R \
         60_ghi_reliability.R \
         61_ghi_scale_and_sensitivity.R \
         63_rob_and_publication_bias.R \
         65_multiplicity.R ; do
  echo "=== $s ==="
  Rscript "$s" || { echo "FAILED: $s"; exit 1; }
done

echo "All scripts completed. Output in $PDI_ROOT/r_analysis/verification/"
