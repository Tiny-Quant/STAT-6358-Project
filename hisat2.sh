#!/bin/bash

QUAL_PHRED=$1
LEN_REQ=$2
TRIM_G=$3
TRIM_X=$4

# Load required modules
module load conda
conda activate rna_seq_env
module load R/4.4.3
module load gcc/11.2.0
module load boost/1.85.0-qggnqmx

# Settings
THREADS=8
HISAT2_INDEX="./indices/HISAT2/genome"  # Path to HISAT2 index
FASTQ_DIR="./fastq_data"
OUTPUT_DIR="hisat2_results"

mkdir -p "$OUTPUT_DIR" "$FASTQ_DIR"
TEMP_DIR=$(mktemp -d)

JOB_ID="${SLURM_ARRAY_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-$RANDOM}"
PARAMS_STR="Q${QUAL_PHRED}_L${LEN_REQ}_G${TRIM_G}_X${TRIM_X}_J${JOB_ID}"
echo "Running HISAT2 pipeline with parameters: $PARAMS_STR"

# Log start time
START_TIME=$(date +%s)

# Process each FASTQ file
for FILE in "$FASTQ_DIR"/*.fastq; do
    BASENAME=$(basename "$FILE" .fastq)
    echo "Processing $BASENAME"

    # Unique prefix for all intermediate files (includes job id)
    PREFIX="$TEMP_DIR/${BASENAME}_${PARAMS_STR}"

    # Preprocess with fastp
    fastp \
        --in1 "$FILE" \
        --qualified_quality_phred "$QUAL_PHRED" \
        --length_required "$LEN_REQ" \
        $( [ "$TRIM_G" -eq 1 ] && echo "--trim_poly_g" ) \
        $( [ "$TRIM_X" -eq 1 ] && echo "--trim_poly_x" ) \
        --out1 "${PREFIX}_trimmed.fastq" \
        --json "${PREFIX}_fastp.json" \
        --html "${PREFIX}_fastp.html"

    # Align with HISAT2 and save SAM
    hisat2 -p "$THREADS" \
        -x "$HISAT2_INDEX" \
        -U "${PREFIX}_trimmed.fastq" \
        -S "$OUTPUT_DIR/${BASENAME}_${PARAMS_STR}.sam"
done

# Log end time
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

# Log elapsed time
echo "Processing completed for $PARAMS_STR. Total time taken: $ELAPSED_TIME seconds."

rm -rf "$TEMP_DIR"
echo "STAR pipeline complete with config: $PARAMS_STR"