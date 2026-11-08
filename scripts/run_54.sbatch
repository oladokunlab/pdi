#!/bin/bash
#SBATCH --job-name=pdi_54
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
set -euo pipefail
module purge
module load GCC/14.3.0 R/4.5.2
export R_LIBS_USER=/scratch/user/$USER/R/library
export PDI_ROOT=/scratch/user/$USER/pdi
export PDI_WORKBOOK=$PDI_ROOT/scripts/GDI_Metadata.xlsx
cd $PDI_ROOT/scripts/R_scripts
Rscript 54_permanova_five_configs.R
