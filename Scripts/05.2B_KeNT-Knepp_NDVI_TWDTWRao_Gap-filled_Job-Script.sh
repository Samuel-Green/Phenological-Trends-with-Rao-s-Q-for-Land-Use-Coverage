#!/bin/bash
#SBATCH --job-name=KeNT-Knepp_NDVI_TWDTWRao
#SBATCH --output=logs/KeNT_%A_%a.out
#SBATCH --error=logs/KeNT_%A_%a.err

#SBATCH --array=1-750
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=7G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd ~/TWDTW_Paper/Scripts
mkdir -p logs

Rscript "05.2A_Knepp_NDVI_TWDTW-RaoQ_SG-Filtered.R"