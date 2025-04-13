library(tximport)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
salmon_dir <- args[1]
output_file <- args[2]

# Locate Salmon output directories
salmon_dirs <- list.dirs(salmon_dir, recursive = FALSE, full.names = TRUE)

# Prepare tx2gene mapping
tx2gene <- read_csv("tx2gene.csv")

# Read Salmon quant files
files <- file.path(salmon_dirs, "quant.sf")
names(files) <- basename(salmon_dirs)

# Import counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                ignoreTxVersion = T, ignoreAfterBar = T)

# Save only the final count matrix
write.csv(txi$counts, file = output_file, row.names = TRUE)

print("Gene count matrix saved.")
