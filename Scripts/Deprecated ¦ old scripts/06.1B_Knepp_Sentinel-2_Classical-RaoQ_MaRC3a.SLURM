#!/bin/bash
#SBATCH --job-name=20Classic_Knepp_PPI
#SBATCH --output=logs/Knepp_2020_PPI_classical_%A_%a.out
#SBATCH --error=logs/Knepp_2020_PPI_classical_%A_%a.err

#SBATCH --array=1-500
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module purge
module load gnu9/9.4.0
module load R/4.1.2

export YEAR=2020   # <-- change per submission

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript 06.1A_Knepp_Estate_Classical-RaoQ_MaRC3a.R