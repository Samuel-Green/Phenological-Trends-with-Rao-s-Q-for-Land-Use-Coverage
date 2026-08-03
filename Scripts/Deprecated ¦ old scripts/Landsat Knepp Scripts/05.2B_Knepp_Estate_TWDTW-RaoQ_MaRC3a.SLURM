#!/bin/bash
#SBATCH --job-name=knepp_NDVI_twdtw
#SBATCH --output=logs/knepp_NDVI_twdtw_%A_%a.out
#SBATCH --error=logs/knepp_NDVI_twdtw_%A_%a.err

#SBATCH --array=1-1000
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G

module purge
module load gnu9/9.4.0
module load R/4.1.2

export YEAR=2000   # <-- change per submission

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript 05.2A_Knepp_Estate_TWDTW-RaoQ_MaRC3a.R