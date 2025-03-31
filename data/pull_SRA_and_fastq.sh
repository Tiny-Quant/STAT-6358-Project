#!/bin/bash

ACCESSION_LIST="SRR_Acc_List.txt"

mkdir -p SRA
mkdir -p fastq_data

while read -r SRR_ID; do

    echo "Processing SRR ID: $SRR_ID"

    if [ ! -f "SRA/$SRR_ID/$SRR_ID.sra" ]; then 
        prefetch $SRR_ID --output-directory SRA
    fi 

    if [ ! -f "fastq_data/$SRR_ID.fastq" ]; then 
        fasterq-dump SRA/$SRR_ID/$SRR_ID.sra --outdir fastq_data --split-files --progress 
    fi

done < "$ACCESSION_LIST"

echo "Processing complete!"
