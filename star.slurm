#!/bin/bash
#SBATCH --job-name=STAR_random
#SBATCH --output=logs/star_random_%A_%a.out
#SBATCH --array=1-10
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --partition=htc
#SBATCH --account=mmcgee_stat6358_1252_001_0001


module load conda 
module load gcc/11.2.0
module load boost/1.85.0-qggnqmx

# Activate the conda environment and load R
conda activate rna_seq_env
module load R/4.4.3


bash star.sh  # internally handles random sampling
