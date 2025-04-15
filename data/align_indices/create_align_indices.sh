#!/bin/bash
set -euo pipefail # Exit on any fail. 

# ----------------------------
# CONFIGURATION
# ----------------------------
THREADS=8
SPECIES="human"
GENCODE_VERSION="44"
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_${SPECIES}/release_${GENCODE_VERSION}"

# Filenames on GENCODE FTP
FASTA_GZ="gencode.v${GENCODE_VERSION}.transcripts.fa.gz"
GTF_GZ="gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
GENOME_FASTA_GZ="GRCh38.primary_assembly.genome.fa.gz"
GENOME_FASTA_URL="${BASE_URL}/${GENOME_FASTA_GZ}"

# ----------------------------
# PATH SETUP
# ----------------------------
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME_DIR="${BASE_DIR}/genome_ref"
SALMON_IDX="${BASE_DIR}/salmon_index"
KALLISTO_IDX="${BASE_DIR}/kallisto_index.idx"
RSEM_REF="${BASE_DIR}/rsem_ref"
TX2GENE_CSV="${BASE_DIR}/tx2gene.csv"

mkdir -p "${GENOME_DIR}" "${SALMON_IDX}" "${BASE_DIR}" "${RSEM_REF}"

echo "=== Reference and index directory: ${BASE_DIR} ==="

# ----------------------------
# 1. Download & unpack GENCODE transcripts
# ----------------------------
cd "${GENOME_DIR}"
echo "-- GENCODE directory: ${GENOME_DIR} --"

if [ ! -f "${FASTA_GZ}" ]; then
  echo "Downloading ${FASTA_GZ}..."
  wget -c "${BASE_URL}/${FASTA_GZ}"
else
  echo "Skipping download of ${FASTA_GZ}, already exists."
fi

if [ ! -f "${GTF_GZ}" ]; then
  echo "Downloading ${GTF_GZ}..."
  wget -c "${BASE_URL}/${GTF_GZ}"
else
  echo "Skipping download of ${GTF_GZ}, already exists."
fi

if [ ! -f "${FASTA_GZ%.gz}" ]; then
  echo "Unpacking ${FASTA_GZ}..."
  gunzip -kf "${FASTA_GZ}"
else
  echo "Skipping unpack of ${FASTA_GZ}, already unpacked."
fi

if [ ! -f "${GTF_GZ%.gz}" ]; then
  echo "Unpacking ${GTF_GZ}..."
  gunzip -kf "${GTF_GZ}"
else
  echo "Skipping unpack of ${GTF_GZ}, already unpacked."
fi

TRANSCRIPTS_FA="${GENOME_DIR}/${FASTA_GZ%.gz}"
ANNOTATION_GTF="${GENOME_DIR}/${GTF_GZ%.gz}"

# ----------------------------
# 2. Build Salmon index
# ----------------------------
if [ ! -f "${SALMON_IDX}/versionInfo.json" ]; then
  echo "Building Salmon index..."
  salmon index \
    --threads ${THREADS} \
    -t "${TRANSCRIPTS_FA}" \
    -i "${SALMON_IDX}"
else
  echo "Skipping Salmon index; already exists at ${SALMON_IDX}."
fi

# ----------------------------
# 3. Build Kallisto index
# ----------------------------
if [ ! -f "${KALLISTO_IDX}" ]; then
  echo "Building Kallisto index..."
  kallisto index \
    -i "${KALLISTO_IDX}" \
    "${TRANSCRIPTS_FA}"
else
  echo "Skipping Kallisto index; already exists at ${KALLISTO_IDX}."
fi

# ----------------------------
# 4. Prepare RSEM (Bowtie2) reference
# ----------------------------
# if [ ! -f "${GENOME_FASTA_GZ}" ]; then
#   echo "Downloading genome FASTA ${GENOME_FASTA_GZ}..."
#   wget -c "${GENOME_FASTA_URL}"
# else
#   echo "Skipping download of ${GENOME_FASTA_GZ}, already exists."
# fi

# if [ ! -f "${GENOME_FASTA_GZ%.gz}" ]; then
#   echo "Unpacking ${GENOME_FASTA_GZ}..."
#   gunzip -kf "${GENOME_FASTA_GZ}"
# else
#   echo "Skipping unpack of ${GENOME_FASTA_GZ}, already unpacked."
# fi

# GENOME_FA="${GENOME_DIR}/${GENOME_FASTA_GZ%.gz}"

# if [ ! -f "${RSEM_REF}/rsem.transcripts.fa" ]; then
#   echo "Preparing RSEM reference (Bowtie2)..."
#   rsem-prepare-reference \
#     --bowtie2 \
#     --gtf "${ANNOTATION_GTF}" \
#     --num-threads ${THREADS} \
#     "${GENOME_FA}" \
#     "${RSEM_REF}/rsem"
# else
#   echo "Skipping RSEM reference; already exists at ${RSEM_REF}."
# fi

# ----------------------------
# 5. Create tx2gene.csv
# ----------------------------
if [ ! -f "$TX2GENE_CSV" ]; then
    echo "Generating tx2gene.csv..."
    gunzip -c "$GTF_GZ" | awk '$3 == "transcript" { 
        split($12, a, "."); 
        split($10, b, "."); 
        print a[1]","b[1]
    }' | tr -d '";' > "$TX2GENE_CSV"
else
    echo "tx2gene.csv already exists. Skipping."
fi

echo "=== All indices built or verified in ${BASE_DIR} ==="
echo "  • Salmon:    ${SALMON_IDX}"
echo "  • Kallisto:  ${KALLISTO_IDX}"
echo "  • RSEM:      ${RSEM_REF}"
echo "  • tx2gene:   ${TX2GENE_CSV}"
