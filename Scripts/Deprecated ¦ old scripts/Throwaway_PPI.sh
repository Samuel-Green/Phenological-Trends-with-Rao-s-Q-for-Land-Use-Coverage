#!/bin/bash
#SBATCH --job-name=PPI_No.4
#SBATCH --partition=normal
#SBATCH --cpus-per-task=40
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output="slurm-Throwaway_PPI-%j.out"

module purge
module load gnu9/9.4.0
module load R/4.1.2

cd "$HOME/TWDTW_Paper/Scripts"

Rscript "Throwaway_Kili_PPI_Calc.R"