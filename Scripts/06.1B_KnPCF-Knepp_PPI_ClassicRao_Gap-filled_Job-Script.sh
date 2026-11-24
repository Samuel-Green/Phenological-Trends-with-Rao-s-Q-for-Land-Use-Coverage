#!/bin/bash
#SBATCH --job-name=KnPCF-Knepp_PPI_ClassicRao
#SBATCH --output=logs/KnPCF_%A_%a.out
#SBATCH --error=logs/KnPCF_%A_%a.err

#SBATCH --array=1-750
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd ~/TWDTW_Paper/Scripts
mkdir -p logs

Rscript "06.1A_Knepp_PPI_Classic-RaoQ_SG-Filtered.R"