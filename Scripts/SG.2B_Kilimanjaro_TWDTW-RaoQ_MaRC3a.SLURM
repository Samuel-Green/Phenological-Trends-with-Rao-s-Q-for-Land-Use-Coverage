#!/bin/bash
#SBATCH --job-name=Kili_SG_TWDTW_RaoQ
#SBATCH --output=logs/SG_tiles/Kili_SG_TWDTW_RaoQ_%A_%a.out
#SBATCH --error=logs/SG_tiles/Kili_SG_TWDTW_RaoQ_%A_%a.err

#SBATCH --array=1-2000

#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript SG.2A_Kilimanjaro_TWDTW-RaoQ_MaRC3a.R