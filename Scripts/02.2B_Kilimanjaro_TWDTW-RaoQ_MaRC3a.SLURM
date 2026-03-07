#!/bin/bash
#SBATCH --job-name=kili_twdtw
#SBATCH --output=logs/kili_%A_%a.out
#SBATCH --error=logs/kili_%A_%a.err

#SBATCH --array=1-72

#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module purge
module load R/4.1.2

Rscript 02.2A_Kilimanjaro_TWDTW-RaoQ_MaRC3a.R