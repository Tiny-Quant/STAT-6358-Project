library(dplyr)
library(readr)
library(tidyr)
library(stringr)

input_dir <- "hisat2_count_matrices"
output_dir <- "hisat2_combined_counts"
if (!dir.exists(output_dir)) dir.create(output_dir)

# List all count files
count_files <- list.files(input_dir, pattern = "_count_matrix.csv$", full.names = TRUE)

# Extract configuration from filename
get_config <- function(filename) {
  str_match(basename(filename), ".*_(Q[0-9]+_L[0-9]+_G[0-9]+_X[0-9]+_J[0-9_]+)_count_matrix.csv")[,2]
}

configs <- unique(sapply(count_files, get_config))

for (config in configs) {
  files <- count_files[sapply(count_files, function(f) get_config(f) == config)]
  count_list <- lapply(files, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    df <- df %>% select(gene_id, counts, sample)
    colnames(df)[2] <- df$sample[1]
    df %>% select(-sample)
  })
  combined <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), count_list)
  combined <- combined %>% arrange(gene_id)
  output_file <- file.path(output_dir, paste0(config, "_matrix.csv"))
  write_csv(combined, output_file)
  cat("✅ Combined HISAT2 matrix saved to:", output_file, "\n")
}
