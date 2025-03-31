#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(ShortRead)
  library(Rsubread)
})

usage <- function(){
  cat("\nUsage: rna_seq_pipeline.R [options]\n\nOptions:\n  --fastq        Path to input FASTQ file (required)\n  --index        Path to reference genome index or index directory (required)\n  --aligner      STAR, HISAT2, Salmon, Kallisto, or Rsubread (default: Rsubread)\n  --output_dir   Directory to save output files (default: './output')\n\n")
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) usage()

params <- list(
  fastq = NULL, index = NULL, aligner = "Rsubread", output_dir = "./output"
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

params$aligner <- tolower(params$aligner)

if(!dir.exists(params$output_dir)){ dir.create(params$output_dir, recursive = TRUE) }

cat("Running pipeline with parameters:\n"); print(params)

### Define Paths to Tools ###
hisat2_bin <- "/projects/dheitjan/hschuler_dissertatio/main_storage/.conda/envs/rna_seq_env/bin/hisat2"
star_bin <- "/projects/dheitjan/hschuler_dissertatio/main_storage/.conda/envs/rna_seq_env/bin/STAR"
salmon_bin <- "/projects/dheitjan/hschuler_dissertatio/main_storage/.conda/envs/rna_seq_env/bin/salmon"
kallisto_bin <- "/projects/dheitjan/hschuler_dissertatio/main_storage/.conda/envs/rna_seq_env/bin/kallisto"
samtools_bin <- "/projects/dheitjan/hschuler_dissertatio/main_storage/.conda/envs/rna_seq_env/bin/samtools"

### ALIGNMENT FUNCTION (FOR STAR, HISAT2, RSUBREAD) ###
run_alignment <- function(fastq_file, index_base, output_dir, aligner){
  if(aligner == "rsubread"){
    bam_output <- file.path(output_dir, paste0(basename(fastq_file), ".bam"))
    cat("Running Rsubread alignment on:", fastq_file, "\n")
    align(index = index_base, readfile1 = fastq_file, output_file = bam_output, nthreads = 4)
    mapped_reads <- countBam(bam_output)$records
  
  } else if(aligner == "star"){
    out_prefix <- file.path(output_dir, paste0(basename(fastq_file), "_STAR_"))
    cmd <- sprintf("%s --runThreadN 4 --genomeDir %s --readFilesIn %s --outFileNamePrefix %s", 
                   star_bin, index_base, fastq_file, out_prefix)
    system(cmd)
    log_file <- file.path(output_dir, paste0(basename(fastq_file), "_STAR_Log.final.out"))
    mapped_reads <- as.numeric(system(sprintf("grep 'Uniquely mapped reads %%' %s | awk '{print $6}'", log_file), intern=TRUE))
  
  } else if(aligner == "hisat2"){
    sam_output <- file.path(output_dir, paste0(basename(fastq_file), "_HISAT2.sam"))
    cmd <- sprintf("%s -p 4 -x %s -U %s -S %s", hisat2_bin, index_base, fastq_file, sam_output)
    system(cmd)
    mapped_reads <- as.numeric(system(sprintf("%s view -c -F 4 %s", samtools_bin, sam_output), intern=TRUE))
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
    cmd <- sprintf("module load gcc/11.2.0; module load boost/1.85.0-qggnqmx; salmon quant -i %s -l A -r %s -p 4 -o %s", 
                   index_base, fastq_file, quant_dir)
  } else if(aligner == "kallisto"){
    cmd <- sprintf("kallisto quant -i %s -o %s %s", index_base, quant_dir, fastq_file)
  }

  system(cmd)

  quant_file <- file.path(quant_dir, "quant.sf")
  if (!file.exists(quant_file)) {
    stop("Error: Quantification output file not found. Check module dependencies and paths.")
  }

  quant_data <- read.table(quant_file, header=TRUE)
  tpm_variance <- var(quant_data$TPM)

  log_entry <- sprintf("%s, %s, %.5f\n", basename(fastq_file), aligner, tpm_variance)
  write(log_entry, file = file.path(output_dir, "expression_log.csv"), append = TRUE)

  unlink(quant_dir, recursive = TRUE)  # Delete entire quantification directory
  return(tpm_variance)
}
### RUN PIPELINE ###
if(params$aligner %in% c("star", "hisat2", "rsubread")){
  mapped_result <- run_alignment(params$fastq, params$index, params$output_dir, params$aligner)
} else if(params$aligner %in% c("salmon", "kallisto")){
  expression_var <- run_quantification(params$fastq, params$index, params$output_dir, params$aligner)
}

cat("Pipeline execution complete.\n")
