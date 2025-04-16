# Load required packages
library(Rsubread)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript HISAT2_count_matrices.R <bam_file> <gtf_file>")
bam_file <- args[1]
gtf_file <- args[2]

output_dir <- "hisat2_count_matrices"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Extract sample ID and configuration from BAM filename
sample_id <- sub("_Q.*", "", basename(bam_file))
configuration <- sub(".*_(Q[0-9]+_L[0-9]+_G[0-9]+_X[0-9]+_J[0-9_]+)\\.bam", "\\1", basename(bam_file))
config_parts <- unlist(strsplit(configuration, "_"))
quality <- as.integer(sub("Q", "", config_parts[1]))
length <- as.integer(sub("L", "", config_parts[2]))
geneFlag <- as.integer(sub("G", "", config_parts[3]))
xFlag <- as.integer(sub("X", "", config_parts[4]))
jobFlag <- sub("J", "", config_parts[5])

count_data <- featureCounts(
  bam_file,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  allowMultiOverlap = TRUE,
  strandSpecific = 0
)

counts <- count_data$counts
count_df <- data.frame(
  gene_id = rownames(counts),
  counts = as.numeric(counts),
  sample = sample_id,
  quality = quality,
  length = length,
  geneFlag = geneFlag,
  xFlag = xFlag,
  jobFlag = jobFlag
)

output_file <- file.path(output_dir, paste0(sample_id, "_", configuration, "_count_matrix.csv"))
write_csv(count_df, output_file)
cat("✅ HISAT2 gene count matrix saved to:", output_file, "\n")