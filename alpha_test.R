#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(ShortRead)
  library(Rsubread)
})

# Function to display usage info
usage <- function(){
  cat("
Usage: rna_seq_pipeline.R [options]

Options:
  --fastq        Path to input FASTQ file (required)
  --threshold    Quality threshold for filtering (default: 20)
  --min_length   Minimum read length after filtering (default: 30)
  --index        Path to reference genome index or index directory (required)
  --aligner      Alignment method: STAR, HISAT2, Salmon, Kallisto, or Rsubread (default: Rsubread)
  --output_dir   Directory to save output files (default: './output')
\n")
  quit(status = 1)
}

# Parse command-line arguments using commandArgs
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) usage()

# Default parameter values
params <- list(
  fastq = NULL,
  threshold = 20,
  min_length = 30,
  index = NULL,
  aligner = "Rsubread",
  output_dir = "./output"
)

# Simple parser for arguments in the form --key=value
for(arg in args){
  if(grepl("^--", arg)){
    split_arg <- strsplit(sub("^--", "", arg), "=")[[1]]
    if(length(split_arg) != 2) usage()
    key <- split_arg[1]
    value <- split_arg[2]
    if(key %in% names(params)){
      params[[key]] <- value
    } else {
      cat("Unknown parameter:", key, "\n")
      usage()
    }
  }
}

# Check for required parameters
if(is.null(params$fastq) || is.null(params$index)){
  cat("Error: --fastq and --index parameters are required.\n")
  usage()
}

# Convert numeric parameters
params$threshold <- as.numeric(params$threshold)
params$min_length <- as.numeric(params$min_length)

# Standardize aligner string to lower-case for comparison
params$aligner <- tolower(params$aligner)

# Create output directory if it doesn't exist
if(!dir.exists(params$output_dir)){
  dir.create(params$output_dir, recursive = TRUE)
}

cat("Running pipeline with parameters:\n")
print(params)

# Define function to trim reads based on quality threshold and minimum length
trim_reads <- function(file, quality_threshold, min_length, output_dir){
  cat("Reading FASTQ file:", file, "\n")
  fq <- readFastq(file)
  scores <- alphabetScore(quality(fq))
  fq_filtered <- fq[scores >= quality_threshold & width(fq) >= min_length]
  
  output_file <- file.path(output_dir, paste0(basename(file), "_trimmed.fastq"))
  writeFastq(fq_filtered, output_file, compress = FALSE)
  cat("Trimmed FASTQ saved to:", output_file, "\n")
  return(output_file)
}

# Define function to run alignment based on selected method
run_alignment <- function(fastq_file, index_base, output_dir, aligner){
  # Prepare an output file prefix or directory name based on aligner
  if(aligner == "rsubread"){
    bam_output <- file.path(output_dir, paste0(basename(fastq_file), ".bam"))
    cat("Running alignment using Rsubread on file:", fastq_file, "\n")
    align(
      index = index_base,
      readfile1 = fastq_file,
      output_file = bam_output,
      nthreads = 4  # Adjust thread count as needed
    )
    cat("Alignment complete. BAM file saved to:", bam_output, "\n")
    return(bam_output)
  } else if(aligner == "star"){
    # Construct STAR command; adjust options and paths as needed
    out_prefix <- file.path(output_dir, paste0(basename(fastq_file), "_STAR_"))
    cmd <- sprintf("STAR --runThreadN 4 --genomeDir %s --readFilesIn %s --outFileNamePrefix %s",
                   index_base, fastq_file, out_prefix)
    cat("Running alignment using STAR:\n", cmd, "\n")
    system(cmd)
    # STAR outputs multiple files; returning the main BAM file path:
    bam_file <- file.path(output_dir, paste0(basename(fastq_file), "_STAR_Aligned.out.bam"))
    cat("STAR alignment complete. BAM file should be at:", bam_file, "\n")
    return(bam_file)
  } else if(aligner == "hisat2"){
    # Construct HISAT2 command; adjust options and paths as needed
    sam_output <- file.path(output_dir, paste0(basename(fastq_file), "_HISAT2.sam"))
    cmd <- sprintf("hisat2 -p 4 -x %s -U %s -S %s",
                   index_base, fastq_file, sam_output)
    cat("Running alignment using HISAT2:\n", cmd, "\n")
    system(cmd)
    # Convert SAM to BAM using samtools (assumed to be installed)
    bam_output <- sub("\\.sam$", ".bam", sam_output)
    convert_cmd <- sprintf("samtools view -Sb %s > %s", sam_output, bam_output)
    system(convert_cmd)
    cat("HISAT2 alignment complete. BAM file saved to:", bam_output, "\n")
    return(bam_output)
  } else if(aligner == "salmon"){
    # Salmon performs quasi-mapping and quantification. Here we run quantification only.
    quant_dir <- file.path(output_dir, paste0(basename(fastq_file), "_Salmon_quant"))
    cmd <- sprintf("salmon quant -i %s -l A -r %s -p 4 -o %s",
                   index_base, fastq_file, quant_dir)
    cat("Running quantification using Salmon:\n", cmd, "\n")
    system(cmd)
    cat("Salmon quantification complete. Results saved to:", quant_dir, "\n")
    return(quant_dir)
  } else if(aligner == "kallisto"){
    # Kallisto quantification command; adjust options as needed
    quant_dir <- file.path(output_dir, paste0(basename(fastq_file), "_Kallisto_quant"))
    cmd <- sprintf("kallisto quant -i %s -o %s %s",
                   index_base, quant_dir, fastq_file)
    cat("Running quantification using Kallisto:\n", cmd, "\n")
    system(cmd)
    cat("Kallisto quantification complete. Results saved to:", quant_dir, "\n")
    return(quant_dir)
  } else {
    stop("Unsupported aligner selected. Please choose from STAR, HISAT2, Salmon, Kallisto, or Rsubread.")
  }
}

# Step 1: Base filtering (trimming)
trimmed_fastq <- trim_reads(params$fastq, params$threshold, params$min_length, params$output_dir)

# Step 2: Alignment / Quantification based on selected aligner
alignment_result <- run_alignment(trimmed_fastq, params$index, params$output_dir, params$aligner)

cat("Pipeline execution complete. Alignment result located at:\n")
print(alignment_result)
