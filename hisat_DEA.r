library(dplyr)
library(tidyr)
library(readr)
library(edgeR)
library(stringr)

set.seed(42) # For reproducibility of random norm_method selection

if (!dir.exists("hisat2_DE_results")) dir.create("hisat2_DE_results")

# Directory containing all HISAT2 count files
count_dir <- "hisat2_count_matrices"
count_files <- list.files(count_dir, pattern = "_count_matrix.csv$", full.names = TRUE)

# Extract pipeline parameters from filenames
extract_params <- function(filename) {
  m <- str_match(basename(filename),
                 "SRR\\d+_Q(\\d+)_L(\\d+)_G(\\d+)_X(\\d+)_J(\\d+)_\\d+_count_matrix.csv"
  )
  data.frame(
    sample = str_extract(basename(filename), "SRR\\d+"),
    quality = m[2],
    length = m[3],
    geneFlag = m[4],
    xFlag = m[5],
    jobFlag = m[6],
    filename = filename,
    stringsAsFactors = FALSE
  )
}

params_df <- bind_rows(lapply(count_files, extract_params))

# Identify all unique pipeline runs
pipeline_runs <- params_df %>%
  select(quality, length, geneFlag, xFlag, jobFlag) %>%
  distinct()

# Treatments vector
treatments <- c(
  "chemo-naive",
  "chemo-naive",
  "chemo-naive",
  "chemo-resistant",
  "chemo-resistant",
  "chemo-resistant",
  "chemo-resistant_prmt5i",
  "chemo-resistant_prmt5i",
  "chemo-resistant_prmt5i"
)

# Define normalization methods to sample from
norm_methods <- c("none", "TMM", "TMMwsp", "RLE", "upperquartile")

# Loop over each pipeline run
for (i in seq_len(nrow(pipeline_runs))) {
  run <- pipeline_runs[i, ]
  run_files <- params_df %>%
    filter(
      quality == run$quality,
      length == run$length,
      geneFlag == run$geneFlag,
      xFlag == run$xFlag,
      jobFlag == run$jobFlag
    )
  
  # Read and merge counts for this run
  count_list <- lapply(run_files$filename, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    df <- df %>% select(gene_id, counts, sample)
    colnames(df)[2] <- df$sample[1] # Name column by sample
    df %>% select(-sample)
  })
  
  # Merge all by gene_id
  count_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), count_list)
  count_matrix <- as.data.frame(count_matrix) # Ensure it's a data.frame, not a tibble
  rownames(count_matrix) <- count_matrix$gene_id
  count_matrix <- count_matrix[, setdiff(colnames(count_matrix), "gene_id")]
  count_matrix <- as.matrix(count_matrix)
  storage.mode(count_matrix) <- "integer"
  
  # Ensure columns are in the same order as treatments
  count_matrix <- count_matrix[, order(colnames(count_matrix))]
  
  # Randomly select a normalization method for this run
  norm_method <- sample(norm_methods, 1)
  
  # Run edgeR
  dge <- DGEList(counts = count_matrix, group = treatments)
  dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
  dge <- estimateDisp(dge)
  if (norm_method != "none") {
    dge <- calcNormFactors(dge, method = norm_method)
  }
  group <- factor(treatments)
  design <- model.matrix(~ group)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef = 2)
  de_results <- topTags(lrt, n = Inf)$table
  
  # Add gene names as a column and write with row.names = FALSE
  de_results$gene <- rownames(de_results)
  de_results <- de_results[, c("gene", setdiff(colnames(de_results), "gene"))]
  
  out_name <- sprintf(
    "DE_hisat2_Q%s_L%s_G%s_X%s_J%s_%s.csv",
    run$quality, run$length, run$geneFlag, run$xFlag, run$jobFlag, norm_method
  )
  write.csv(de_results, file = file.path("hisat2_DE_results", out_name), row.names = FALSE)
  
  cat("Finished pipeline run:", out_name, "\n")
}