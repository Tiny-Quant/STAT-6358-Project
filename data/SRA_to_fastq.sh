#!/bin/bash

ACCESSION_LIST="SRR_Acc_List.txt"

mkdir -p SRA
mkdir -p dump_1 dump_2 dump_3

while read -r SRR_ID; do

    echo "Processing SRR ID: $SRR_ID"

    prefetch $SRR_ID --output-directory SRA

    # Option 1: Default
    fasterq-dump SRA/$SRR_ID/$SRR_ID.sra --outdir dump_1 --split-files --progress 

    # Option 2: Skip Technical 
    fasterq-dump SRA/$SRR_ID/$SRR_ID.sra --outdir dump_2 --split-files --progress --skip-technical

    # Option 3: 10 Threads
    fasterq-dump SRA/$SRR_ID/$SRR_ID.sra --outdir dump_3 --split-files --progress --threads 10

done < "$ACCESSION_LIST"

echo "Processing complete!"

