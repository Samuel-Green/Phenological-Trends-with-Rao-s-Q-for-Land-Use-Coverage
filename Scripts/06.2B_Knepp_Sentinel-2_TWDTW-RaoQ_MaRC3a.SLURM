#!/bin/bash
#SBATCH --job-name=20twdtw_Knepp_PPI
#SBATCH --output=logs/knepp_2020_PPI_twdtw_%A_%a.out
#SBATCH --error=logs/knepp_2020_PPI_twdtw_%A_%a.err

#SBATCH --array=1-500
#SBATCH --time=00:55:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G

module purge
module load gnu9/9.4.0
module load R/4.1.2

export YEAR=2020   # <-- change per submission

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript 06.2A_Knepp_Estate_TWDTW-RaoQ_MaRC3a.R