#!/bin/bash

THREADS=8
SALMON_INDEX="salmon_index"
TX2GENE="tx2gene.csv"
SRR_LIST="../SRR_Acc_List.txt"  # File containing SRR IDs, one per line
FASTQ_DIR="../fastq_data"
OUTPUT_DIR="salmon_count_matrices"
R_SCRIPT="generate_salmon_count_matrix.R"

mkdir -p "$OUTPUT_DIR" "$FASTQ_DIR"

TEMP_DIR=$(mktemp -d)

PARAMS_STR=""
# Randomly select parameters. 
QUAL_PHRED=$(awk -v min=20 -v max=30 'BEGIN{srand(); print int(min+rand()*(max-min+1))}')
LEN_REQ=$(awk -v min=30 -v max=50 'BEGIN{srand(); print int(min+rand()*(max-min+1))}')

TRIM_G=$((RANDOM % 2))
TRIM_X=$((RANDOM % 2))

PARAMS_STR="Q${QUAL_PHRED}_L${LEN_REQ}_G${TRIM_G}_X${TRIM_X}"

echo "Processed files with QUAL_PHRED=$QUAL_PHRED, LEN_REQ=$LEN_REQ, TRIM_G=$TRIM_G, TRIM_X=$TRIM_X"

conda install -c bioconda fastp salmon

for FILE in "$FASTQ_DIR"/*.fastq; do 

    echo "Now processing $FILE..."

    BASENAME=$(basename "$FILE" .fastq)

    # Filter and trim based on sampled parameters. 
    fastp \
        --in1 "$FILE" \
        --qualified_quality_phred "$QUAL_PHRED" \
        --length_required "$LEN_REQ" \
        $( [ "$TRIM_G" -eq 1 ] && echo "--trim_poly_g" ) \
        $( [ "$TRIM_X" -eq 1 ] && echo "--trim_poly_x" ) \
        --out1 "$TEMP_DIR/${BASENAME}_trimmed.fastq" \
        --json "$TEMP_DIR/${BASENAME}_fastp.json" \
        --html "$TEMP_DIR/${BASENAME}_fastp.html"

    echo "hello" 

    # Run Salmon quantification
    salmon quant -i "$SALMON_INDEX" \
        -l A \
        -r "$TEMP_DIR/${BASENAME}_trimmed.fastq" \
        -o "$TEMP_DIR/${BASENAME}_salmon" \
        --gcBias --seqBias --validateMappings

done

# Define output file with selected parameters
FINAL_COUNT_MATRIX="$OUTPUT_DIR/gene_count_matrix_${PARAMS_STR}.csv"

# Run R script to generate gene count matrix
Rscript "$R_SCRIPT" "$TEMP_DIR" "$FINAL_COUNT_MATRIX"

# Remove temporary files
rm -rf "$TEMP_DIR"

echo "Pipeline complete. Final count matrix stored in $FINAL_COUNT_MATRIX"