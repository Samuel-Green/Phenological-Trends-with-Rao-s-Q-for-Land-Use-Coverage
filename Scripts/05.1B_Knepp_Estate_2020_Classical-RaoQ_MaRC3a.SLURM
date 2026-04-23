#!/bin/bash
#SBATCH --job-name=2020_classical_knepp_NDVI
#SBATCH --output=logs/knepp_2020_NDVI_classical_%A_%a.out
#SBATCH --error=logs/knepp_2020_NDVI_classical_%A_%a.err

#SBATCH --array=1-500
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module purge
module load gnu9/9.4.0
module load R/4.1.2

export YEAR=2020   # <-- change per submission

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript 05.1A_Knepp_Estate_Classical-RaoQ_MaRC3a.R