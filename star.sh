#!/bin/bash

# Load required modules
module load conda
conda activate rna_seq_env
module load R/4.4.3
module load gcc/11.2.0
module load boost/1.85.0-qggnqmx

THREADS=8
STAR_INDEX="./indices/STAR"  # update path
FASTQ_DIR="./fastq_data"
OUTPUT_DIR="star_results"


mkdir -p "$OUTPUT_DIR" "$FASTQ_DIR"
TEMP_DIR=$(mktemp -d)

# Randomly select QC parameters
QUAL_PHRED=$(awk -v min=20 -v max=30 'BEGIN{srand(); print int(min+rand()*(max-min+1))}')
LEN_REQ=$(awk -v min=30 -v max=50 'BEGIN{srand(); print int(min+rand()*(max-min+1))}')
TRIM_G=$((RANDOM % 2))
TRIM_X=$((RANDOM % 2))

# Get SLURM job and task ID or fallback to process ID for uniqueness
JOB_ID="${SLURM_ARRAY_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-$RANDOM}"

PARAMS_STR="Q${QUAL_PHRED}_L${LEN_REQ}_G${TRIM_G}_X${TRIM_X}_J${JOB_ID}"
echo "Running STAR pipeline with: $PARAMS_STR"

for FILE in "$FASTQ_DIR"/*.fastq; do
    BASENAME=$(basename "$FILE" .fastq)
    echo "Processing $BASENAME"

    # Unique prefix for all intermediate files (includes job id)
    PREFIX="$TEMP_DIR/${BASENAME}_${PARAMS_STR}"

    # Trim & filter
    fastp \
        --in1 "$FILE" \
        --qualified_quality_phred "$QUAL_PHRED" \
        --length_required "$LEN_REQ" \
        $( [ "$TRIM_G" -eq 1 ] && echo "--trim_poly_g" ) \
        $( [ "$TRIM_X" -eq 1 ] && echo "--trim_poly_x" ) \
        --out1 "${PREFIX}_trimmed.fastq" \
        --json "${PREFIX}_fastp.json" \
        --html "${PREFIX}_fastp.html"

    # Run STAR
    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "${PREFIX}_trimmed.fastq" \
        --outFileNamePrefix "${PREFIX}_" \
        --outSAMtype BAM Unsorted

    # Move BAM to output
    mv "${PREFIX}_Aligned.out.bam" "$OUTPUT_DIR/${BASENAME}_${PARAMS_STR}.bam"
done

rm -rf "$TEMP_DIR"
echo "STAR pipeline complete with config: $PARAMS_STR"