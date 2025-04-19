#!/bin/bash
#SBATCH --job-name=Dual_Aligner_Run
#SBATCH --output=logs/dual_run_%A_%a.out
#SBATCH --array=1-200
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --partition=standard-s
#SBATCH --account=mmcgee_stat6358_1252_001_0001

module load conda 
module load gcc/11.2.0
module load boost/1.85.0-qggnqmx
conda activate rna_seq_env
module load R/4.4.3

LOG_FILE="hs_parameter_log.txt"

# Determine the 0-based index into the log file
if [ "$SLURM_ARRAY_TASK_ID" -le 100 ]; then
    COMBO_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
    ALIGNER="STAR"
else
    COMBO_INDEX=$((SLURM_ARRAY_TASK_ID - 101))
    ALIGNER="HISAT2"
fi

# Extract the correct combo
COMBO=$(sed -n "$((COMBO_INDEX + 1))p" "$LOG_FILE")
IFS="_" read -r QUAL_PHRED LEN_REQ TRIM_G TRIM_X <<< "$COMBO"

echo "?? Running $ALIGNER with combo: $COMBO"

if [ "$ALIGNER" = "STAR" ]; then
    bash star.sh "$QUAL_PHRED" "$LEN_REQ" "$TRIM_G" "$TRIM_X"
else
    bash hisat2.sh "$QUAL_PHRED" "$LEN_REQ" "$TRIM_G" "$TRIM_X"
fi
