#!/bin/bash
set -euo pipefail # Exit on any fail. 

# Configuration
THREADS=9
PHRED_MIN=20
PHRED_MAX=30
LENGTH_MIN=30
LENGTH_MAX=50
#ALIGNERS=(salmon kallisto rsem)
ALIGNERS=(salmon kallisto)

# File paths. 
SALMON_INDEX="../align_indices/salmon_index"
KALLISTO_INDEX="../align_indices/kallisto_index.idx"
#RSEM_INDEX="../align_indices/rsem_ref"
FASTQ_DIR="../fastq_data"

#TEMP_DIR=$(mktemp -d)
#TEMP_DIR="test_dir"
#mkdir -p "$TEMP_DIR"
TEMP_DIR="/lustre/scratch/client/users/ataychameekiatchai/tmp_parallel_$$"
mkdir -p "$TEMP_DIR"
trap 'echo "Cleaning up $TEMP_DIR"; rm -rf "$TEMP_DIR"' EXIT # Clear files on any script exit. 

# Sample pipeline parameters. 
min_phred=$(( RANDOM % (PHRED_MAX - PHRED_MIN + 1) + PHRED_MIN ))
min_length=$(( RANDOM % (LENGTH_MAX - LENGTH_MIN + 1) + LENGTH_MIN ))
trim_poly_g=$(( RANDOM % 2 ))
trim_poly_x=$(( RANDOM % 2 ))
#aligner=${ALIGNERS[$(( RANDOM % ${#ALIGNERS[@]} ))]}
aligner="salmon"

# Track start time
start_time=$(date +%s)

process_sample() {
    local FILE=$1
    local BASENAME=$(basename "$FILE" .fastq)

    echo "Now processing $FILE…"

    # 1) trim/filter
    fastp \
      --in1 "$FILE" \
      --qualified_quality_phred "$min_phred" \
      --length_required "$min_length" \
      $( [ "$trim_poly_g" -eq 1 ] && echo "--trim_poly_g" ) \
      $( [ "$trim_poly_x" -eq 1 ] && echo "--trim_poly_x" ) \
      --out1 "$TEMP_DIR/${BASENAME}_trimmed.fastq" \
      --json "$TEMP_DIR/${BASENAME}_fastp.json" \
      --html "$TEMP_DIR/${BASENAME}_fastp.html"

    # 2) align
    case $aligner in
      salmon)
        salmon quant -i "$SALMON_INDEX" \
          -l A \
          -r "$TEMP_DIR/${BASENAME}_trimmed.fastq" \
          -o "$TEMP_DIR/${BASENAME}_salmon" \
          --gcBias --seqBias --validateMappings
        ;;
      kallisto)
        kallisto quant -i "$KALLISTO_INDEX" \
          -o "$TEMP_DIR/${BASENAME}_kallisto" \
          -t $THREADS \
          --single -l 200 -s 20 \
          "$TEMP_DIR/${BASENAME}_trimmed.fastq"
        ;;
    #   rsem)
    #     rsem-calculate-expression \
    #       --bowtie2 \
    #       --num-threads "$THREADS" \
    #       "$TEMP_DIR/${BASENAME}_trimmed.fastq" \
    #       "$RSEM_INDEX/rsem" \
    #       "$TEMP_DIR/${BASENAME}_rsem"
    #     ;;
      *)
        echo "Invalid aligner: $aligner"
        ;;
    esac
}

export -f process_sample
export THREADS PHRED_MIN PHRED_MAX LENGTH_MIN LENGTH_MAX \
       SALMON_INDEX KALLISTO_INDEX RSEM_INDEX TEMP_DIR \
       min_phred min_length trim_poly_g trim_poly_x aligner
export TMPDIR="$TEMP_DIR"

find "$FASTQ_DIR" -name "*.fastq" -print0 | \
    parallel --no-notice -0 -j "$THREADS" \
        --tmpdir "$TEMP_DIR" \
        --compress \
        --will-cite \
        --bar --eta --linebuffer \
        process_sample {}

# Record elapsed time
end_time=$(date +%s)
elapsed_sec=$(( end_time - start_time ))

#export PATH="/usr/local/bin:$PATH"
/usr/local/bin/Rscript "./transform_sample.r" \
    "$TEMP_DIR" \
    "$aligner" \
    "$min_phred" \
    "$min_length" \
    "$trim_poly_g" \
    "$trim_poly_x" \
    "$elapsed_sec"

echo "Pipeline sample complete." 

