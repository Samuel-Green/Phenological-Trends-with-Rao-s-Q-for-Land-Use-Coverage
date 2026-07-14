#!/bin/bash
#SBATCH --job-name=Kili_PPI_ClassicRaoQ
#SBATCH --output=logs/Kili_PPI_ClassicRaoQ_%A_%a.out
#SBATCH --error=logs/Kili_PPI_ClassicRaoQ_%A_%a.err

#SBATCH --array=1-2000

#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript 04.1A_Kilimanjaro_Classical-RaoQ_MaRC3a.R