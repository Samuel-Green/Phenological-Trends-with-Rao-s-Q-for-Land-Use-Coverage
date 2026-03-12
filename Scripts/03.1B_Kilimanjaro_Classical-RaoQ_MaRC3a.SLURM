#!/bin/bash
#SBATCH --job-name=kili_kleine_cube
#SBATCH --output=logs/tiny_tiles/kili_MultiCore_%A_%a.out
#SBATCH --error=logs/tiny_tiles/kili_MultiCore_%A_%a.err

#SBATCH --array=1-2000

#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript 03.1C_Kilimanjaro_Classical-RaoQ_MaRC3a.R