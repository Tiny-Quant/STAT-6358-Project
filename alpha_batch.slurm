#!/bin/bash
#SBATCH --job-name=RNAseq_matrix
#SBATCH --output=logs/alpha_test_%A_%a.out  # Save output logs, one per job
#SBATCH --account=mmcgee_stat6358_1252_001_0001
#SBATCH --time=2:00:00           # Adjusted time allocation
#SBATCH --cpus-per-task=16       # Adjusted CPU allocation
#SBATCH --mem=128G               # Increased memory allocation
#SBATCH --array=1-15             # 3 thresholds x 5 aligners = 15 jobs

module load conda 
module load gcc/11.2.0
module load boost/1.85.0-qggnqmx

# Activate the conda environment and load R
conda activate rna_seq_env
module load R/4.4.3

# Debugging: Print environment details
echo "Current Environment: $(conda info --envs | grep \*)"
echo "PATH: $PATH"
echo "STAR location: $(which STAR)"
echo "HISAT2 location: $(which hisat2)"
echo "Salmon location: $(which salmon)"
echo "Kallisto location: $(which kallisto)"


# Define input FASTQ file and output directory
FASTQ="/lustre/work/client/projects/mmcgee/stat6358_1252_001/practice_folder/OvarianCancer/SRR31476642.fastq"
OUTPUT_DIR="/users/hschuler/art_hunter/STAT-6358-Project/alpha_output"

# Define threshold settings (Phred score) and aligners
THRESHOLDS=(25 30 39)
ALIGNERS=("STAR" "HISAT2" "Salmon" "Kallisto" "Rsubread")

# Define separate index paths for each aligner (adjust as needed)
declare -A INDEX_PATHS
INDEX_PATHS["STAR"]="/users/hschuler/art_hunter/STAT-6358-Project/indices/STAR/"
INDEX_PATHS["HISAT2"]="/users/hschuler/art_hunter/STAT-6358-Project/indices/HISAT2/genome"
INDEX_PATHS["Salmon"]="/users/hschuler/art_hunter/STAT-6358-Project/indices/Salmon/"
INDEX_PATHS["Kallisto"]="/users/hschuler/art_hunter/STAT-6358-Project/indices/Kallisto/transcripts.idx"
INDEX_PATHS["Rsubread"]="/users/hschuler/art_hunter/STAT-6358-Project/indices/Rsubread/genome"

# Compute which threshold and aligner to use based on SLURM_ARRAY_TASK_ID
THRESHOLD_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1) / 5 ))
ALIGNER_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1) % 5 ))
THRESHOLD=${THRESHOLDS[$THRESHOLD_INDEX]}
ALIGNER=${ALIGNERS[$ALIGNER_INDEX]}
INDEX=${INDEX_PATHS[$ALIGNER]}

# Print job settings
echo "Running job $SLURM_ARRAY_TASK_ID with:"
echo "  FASTQ = $FASTQ"
echo "  Quality Threshold = $THRESHOLD"
echo "  Aligner = $ALIGNER"
echo "  Index = $INDEX"

# Create output directory for this run
JOB_OUTPUT_DIR="${OUTPUT_DIR}/thr${THRESHOLD}_${ALIGNER}"
mkdir -p "$JOB_OUTPUT_DIR"

# Run fastq-filter externally to filter the FASTQ file.
# The -q option sets the minimum mean quality (Phred score). For example, -q 30 means Q30.
FILTERED_FASTQ=$(mktemp --suffix=.fastq)
echo "Filtering FASTQ with threshold $THRESHOLD using fastq-filter"
cmd="fastq-filter -q ${THRESHOLD} -o ${FILTERED_FASTQ} ${FASTQ}"
echo "Running command: $cmd"
$cmd

# Compute read counts for logging (assuming FASTQ header lines start with '@')
TOTAL_READS=$(grep -c "^@" "$FASTQ")
FILTERED_READS=$(grep -c "^@" "$FILTERED_FASTQ")
echo "Total reads: $TOTAL_READS, Filtered reads: $FILTERED_READS" >> "$JOB_OUTPUT_DIR/quality_log.txt"

# Run the R pipeline (alignment/quantification) using the filtered FASTQ file.
# Note: No threshold is passed anymore to the R script
Rscript alpha_test.R \
  --fastq="$FILTERED_FASTQ" \
  --index="$INDEX" \
  --aligner="$ALIGNER" \
  --output_dir="$JOB_OUTPUT_DIR"

# Remove the filtered FASTQ to conserve space
rm "$FILTERED_FASTQ"

echo "Job $SLURM_ARRAY_TASK_ID complete!"
