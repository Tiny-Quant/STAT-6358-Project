#!/bin/bash
#SBATCH --job-name=STAR_count
#SBATCH --output=logs/star_count_%A_%a.out
#SBATCH --array=1-80
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --partition=htc
#SBATCH --account=mmcgee_stat6358_1252_001_0001

module load conda
conda activate rna_seq_env
module load R/4.4.3

bam_dir="star_results"
gtf_file="Homo_sapiens.GRCh38.110.gtf"

bam_files=($bam_dir/*.bam)
if [ ${#bam_files[@]} -eq 0 ]; then
  echo "Error: No BAM files found in directory $bam_dir." >&2
  exit 1
fi
if [ -z "$SLURM_ARRAY_TASK_ID" ] || [ "$SLURM_ARRAY_TASK_ID" -lt 1 ] || [ "$SLURM_ARRAY_TASK_ID" -gt "${#bam_files[@]}" ]; then
  echo "Error: SLURM_ARRAY_TASK_ID is out of range or undefined." >&2
  exit 1
fi
bam_file=${bam_files[$SLURM_ARRAY_TASK_ID - 1]}

Rscript --vanilla STAR_count_matrices.R "$bam_file" "$gtf_file"