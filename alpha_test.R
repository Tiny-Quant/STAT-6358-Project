#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(ShortRead)
  library(Rsubread)
})

usage <- function(){
  cat("
Usage: rna_seq_pipeline.R [options]

Options:
  --fastq        Path to input FASTQ file (required)
  --threshold    Quality threshold for filtering (default: 20)
  --min_length   Minimum read length after filtering (default: 30)
  --index        Path to reference genome index or index directory (required)
  --aligner      STAR, HISAT2, Salmon, Kallisto, or Rsubread (default: Rsubread)
  --output_dir   Directory to save output files (default: './output')
\n")
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) usage()

params <- list(
  fastq = NULL, threshold = 20, min_length = 30,
  index = NULL, aligner = "Rsubread", output_dir = "./output"
)

for(arg in args){
  if(grepl("^--", arg)){
    split_arg <- strsplit(sub("^--", "", arg), "=")[[1]]
    if(length(split_arg) != 2) usage()
    key <- split_arg[1]; value <- split_arg[2]
    if(key %in% names(params)){ params[[key]] <- value } else {
      cat("Unknown parameter:", key, "\n"); usage()
    }
  }
}

if(is.null(params$fastq) || is.null(params$index)){
  cat("Error: --fastq and --index are required.\n"); usage()
}

params$threshold <- as.numeric(params$threshold)
params$min_length <- as.numeric(params$min_length)
params$aligner <- tolower(params$aligner)

if(!dir.exists(params$output_dir)){ dir.create(params$output_dir, recursive = TRUE) }

cat("Running pipeline with parameters:\n"); print(params)

### TRIM READS ###
trim_reads <- function(file, quality_threshold, min_length, output_dir){
  cat("Reading FASTQ file:", file, "\n")
  fq <- readFastq(file)
  total_reads <- length(fq)
  scores <- alphabetScore(quality(fq))
  mean_quality <- mean(scores)
  fq_filtered <- fq[scores >= quality_threshold & width(fq) >= min_length]
  filtered_reads <- length(fq_filtered)

  log_entry <- sprintf("%s, %d, %d, %d, %d, %.2f\n",
                       basename(file), quality_threshold, min_length,
                       total_reads, filtered_reads, mean_quality)
  write(log_entry, file = file.path(output_dir, "quality_log.csv"), append = TRUE)

  cat("Filtered reads:", filtered_reads, "\n")
  return(filtered_reads)
}

### ALIGNMENT FUNCTION (FOR STAR, HISAT2, RSUBREAD) ###
run_alignment <- function(fastq_file, index_base, output_dir, aligner){
  if(aligner == "rsubread"){
    bam_output <- file.path(output_dir, paste0(basename(fastq_file), ".bam"))
    cat("Running Rsubread alignment on:", fastq_file, "\n")
    align(index = index_base, readfile1 = fastq_file, output_file = bam_output, nthreads = 4)
    mapped_reads <- countBam(bam_output)$records

  } else if(aligner == "star"){
    out_prefix <- file.path(output_dir, paste0(basename(fastq_file), "_STAR_"))
    cmd <- sprintf("STAR --runThreadN 4 --genomeDir %s --readFilesIn %s --outFileNamePrefix %s",
                   index_base, fastq_file, out_prefix)
    system(cmd)
    log_file <- file.path(output_dir, paste0(basename(fastq_file), "_STAR_Log.final.out"))
    mapped_reads <- as.numeric(system(sprintf("grep 'Uniquely mapped reads %%' %s | awk '{print $6}'", log_file), intern=TRUE))

  } else if(aligner == "hisat2"){
    sam_output <- file.path(output_dir, paste0(basename(fastq_file), "_HISAT2.sam"))
    cmd <- sprintf("hisat2 -p 4 -x %s -U %s -S %s", index_base, fastq_file, sam_output)
    system(cmd)
    mapped_reads <- as.numeric(system(sprintf("samtools view -c -F 4 %s", sam_output), intern=TRUE))
  }

  log_entry <- sprintf("%s, %s, %d\n", basename(fastq_file), aligner, mapped_reads)
  write(log_entry, file = file.path(output_dir, "alignment_log.csv"), append = TRUE)
  
  # Remove large alignment files to save space
  if(aligner == "rsubread"){ file.remove(bam_output) }
  if(aligner == "star"){ file.remove(file.path(output_dir, paste0(basename(fastq_file), "_STAR_Aligned.out.bam"))) }
  if(aligner == "hisat2"){ file.remove(sam_output) }

  return(mapped_reads)
}

### QUANTIFICATION FUNCTION (FOR SALMON & KALLISTO) ###
run_quantification <- function(fastq_file, index_base, output_dir, aligner){
  quant_dir <- file.path(output_dir, paste0(basename(fastq_file), "_", aligner, "_quant"))
  
  if(aligner == "salmon"){
    cmd <- sprintf("salmon quant -i %s -l A -r %s -p 4 -o %s", index_base, fastq_file, quant_dir)
  } else if(aligner == "kallisto"){
    cmd <- sprintf("kallisto quant -i %s -o %s %s", index_base, quant_dir, fastq_file)
  }
  
  system(cmd)
  
  quant_file <- file.path(quant_dir, "quant.sf")
  quant_data <- read.table(quant_file, header=TRUE)
  tpm_variance <- var(quant_data$TPM)

  log_entry <- sprintf("%s, %s, %.5f\n", basename(fastq_file), aligner, tpm_variance)
  write(log_entry, file = file.path(output_dir, "expression_log.csv"), append = TRUE)

  unlink(quant_dir, recursive = TRUE)  # Delete entire quantification directory
  return(tpm_variance)
}

### RUN PIPELINE ###
filtered_reads <- trim_reads(params$fastq, params$threshold, params$min_length, params$output_dir)

if(params$aligner %in% c("star", "hisat2", "rsubread")){
  mapped_result <- run_alignment(params$fastq, params$index, params$output_dir, params$aligner)
} else if(params$aligner %in% c("salmon", "kallisto")){
  expression_var <- run_quantification(params$fastq, params$index, params$output_dir, params$aligner)
}

cat("Pipeline execution complete.\n")
