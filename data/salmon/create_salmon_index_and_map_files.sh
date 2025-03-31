#!/bin/bash

# Set variables
GENCODE_VERSION="44"
GENCODE_BASE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}"
TRANSCRIPTOME="gencode.v${GENCODE_VERSION}.transcripts.fa.gz"
GTF_FILE="gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
TX2GENE="tx2gene.csv"
SALMON_INDEX="salmon_index"

# Download salmon using conda. 
conda install -c bioconda salmon

# Download transcriptome FASTA if missing
if [ ! -f "$TRANSCRIPTOME" ]; then
    echo "Downloading GENCODE transcriptome..."
    wget "${GENCODE_BASE}/${TRANSCRIPTOME}"
fi

# Download GTF annotation file if missing
if [ ! -f "$GTF_FILE" ]; then
    echo "Downloading GENCODE GTF annotation..."
    wget "${GENCODE_BASE}/${GTF_FILE}"
fi

# Build Salmon index if not exists
if [ ! -d "$SALMON_INDEX" ]; then
    echo "Building Salmon index..."
    gunzip -c "$TRANSCRIPTOME" > gencode_transcripts.fa
    salmon index -t gencode_transcripts.fa -i "$SALMON_INDEX" -p 8
    rm gencode_transcripts.fa  # Clean up extracted FASTA
else
    echo "Salmon index already exists. Skipping."
fi

# Generate tx2gene.csv if not exists
if [ ! -f "$TX2GENE" ]; then
    echo "Generating tx2gene.csv..."
    gunzip -c "$GTF_FILE" | awk '$3 == "transcript" { 
        split($12, a, "."); 
        split($10, b, "."); 
        print a[1]","b[1]
    }' | tr -d '";' > "$TX2GENE"
else
    echo "tx2gene.csv already exists. Skipping."
fi

