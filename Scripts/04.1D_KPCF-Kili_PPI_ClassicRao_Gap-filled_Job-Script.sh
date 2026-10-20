#!/bin/bash
#SBATCH --job-name=KPCF-LongRun-Kili_PPI_ClassicRao_Gap-filled
#SBATCH --output=logs/KPCF_%A_%a.out
#SBATCH --error=logs/KPCF_%A_%a.err

#SBATCH --array=1-2000

#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=9G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd ~/TWDTW_Paper/Scripts

Rscript "04.1B_Kilimanjaro_PPI_Classic-RaoQ_SG-Filtered.R"