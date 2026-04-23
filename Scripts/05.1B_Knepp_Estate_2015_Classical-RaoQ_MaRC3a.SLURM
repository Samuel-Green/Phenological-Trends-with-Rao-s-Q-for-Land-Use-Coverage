#!/bin/bash
#SBATCH --job-name=2015_classical_knepp_NDVI # <-- change year with submission
#SBATCH --output=logs/knepp_2015_NDVI_classical_%A_%a.out # <-- change year with submission
#SBATCH --error=logs/knepp_2015_NDVI_classical_%A_%a.err # <-- change year with submission

#SBATCH --array=1-500
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module purge
module load gnu9/9.4.0
module load R/4.1.2

export YEAR=2015   # <-- change per submission

cd "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Scripts"

Rscript 05.1A_Knepp_Estate_Classical-RaoQ_MaRC3a.R