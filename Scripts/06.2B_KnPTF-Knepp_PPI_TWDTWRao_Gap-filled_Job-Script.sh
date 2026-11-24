#!/bin/bash
#SBATCH --job-name=KnPTF-Knepp_PPI_TWDTWRao
#SBATCH --output=logs/KnPTF_%A_%a.out
#SBATCH --error=logs/KnPTF_%A_%a.err

#SBATCH --array=1-750
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd ~/TWDTW_Paper/Scripts
mkdir -p logs

Rscript "06.2A_Knepp_PPI_TWDTW-RaoQ_SG-Filtered.R"