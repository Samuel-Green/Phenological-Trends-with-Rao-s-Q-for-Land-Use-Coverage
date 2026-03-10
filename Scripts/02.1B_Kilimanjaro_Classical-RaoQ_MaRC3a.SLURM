#!/bin/bash
#SBATCH --job-name=kili_rao
#SBATCH --output=logs/kili_rao_%A_%a.out
#SBATCH --error=logs/kili_rao_%A_%a.err

#SBATCH --array=1-72

#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript 02.1A_Kilimanjaro_Classical-RaoQ_MaRC3a.R