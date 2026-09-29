#!/bin/bash
#SBATCH --job-name=KNTF-Kili_NDVI_TWDTWRao_Gap-filled
#SBATCH --output=logs/KNTF_%A_%a.out
#SBATCH --error=logs/KNTF_%A_%a.err

#SBATCH --array=1-2000

#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd ~/TWDTW_Paper/Scripts

Rscript "03.2B_Kilimanjaro_TWDTW-RaoQ_SG-Filtered.R"