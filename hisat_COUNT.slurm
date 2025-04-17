#!/bin/bash
#SBATCH --job-name=HISAT2_count
#SBATCH --output=logs/hisat2_count_%A_%a.out
#SBATCH --array=1-90
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1:00:00
#SBATCH --partition=htc
#SBATCH --account=mmcgee_stat6358_1252_001_0001

module load conda
conda activate rna_seq_env
module load R/4.4.3

sam_dir="hisat2_results"
gtf_file="Homo_sapiens.GRCh38.110.gtf"

sam_files=($sam_dir/*.sam)
if [ ${#sam_files[@]} -eq 0 ]; then
  echo "Error: No SAM files found in directory $sam_dir." >&2
  exit 1
fi
if [ -z "$SLURM_ARRAY_TASK_ID" ] || [ "$SLURM_ARRAY_TASK_ID" -lt 1 ] || [ "$SLURM_ARRAY_TASK_ID" -gt "${#sam_files[@]}" ]; then
  echo "Error: SLURM_ARRAY_TASK_ID is out of range or undefined." >&2
  exit 1
fi
sam_file=${sam_files[$SLURM_ARRAY_TASK_ID - 1]}

bam_file="${sam_file%.sam}.bam"
if [ ! -f "$bam_file" ]; then
  samtools view -bS "$sam_file" > "$bam_file"
fi

Rscript --vanilla HISAT2_count_matrices.R "$bam_file" "$gtf_file"