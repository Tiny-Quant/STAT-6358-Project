library(readr)
library(dplyr)
library(stringr)

# ---- Paths ----
star_count_dir <- "star_combined_counts"
hisat2_count_dir <- "hisat2_combined_counts"
star_de_dir <- "star_DE_results"
hisat2_de_dir <- "hisat2_DE_results"

nih_counts <- read_csv("/lustre/work/client/projects/dheitjan/hschuler_dissertatio/main_storage/STAT-6358-Project/data/gen_samples/nih_counts_wide.csv")
nih_de <- read_csv("/lustre/work/client/projects/dheitjan/hschuler_dissertatio/main_storage/STAT-6358-Project/data/gen_samples/nih_DE.csv")

# ---- Helpers ----

# For DE files: DE_hisat2_Q20_L35_G0_X0_J19069379.csv or DE_star_Q20_L35_G0_X0_J19069379.csv
extract_config_de <- function(fname) {
  m <- str_match(fname, "Q([0-9]+)_L([0-9]+)_G([0-9]+)_X([0-9]+)_J([0-9]+)_([A-Za-z0-9]+)\\.csv$")
  list(
    min_phred = as.integer(m[2]),
    min_length = as.integer(m[3]),
    trim_poly_g = as.integer(m[4]),
    trim_poly_x = as.integer(m[5]),
    jobFlag = m[6],
    norm_method = m[7]
  )
}

# For count files: Q20_L35_G0_X0_J19069379_158_matrix.csv
extract_config_count <- function(fname) {
  m <- str_match(fname, "Q([0-9]+)_L([0-9]+)_G([0-9]+)_X([0-9]+)_J[0-9]+")
  list(
    min_phred = as.integer(m[2]),
    min_length = as.integer(m[3]),
    trim_poly_g = as.integer(m[4]),
    trim_poly_x = as.integer(m[5])
  )
}

summarize_counts <- function(file, nih_counts) {
  df <- read_csv(file, show_col_types = FALSE)
  colnames(df)[1] <- "gene_id"
  common_genes <- intersect(df$gene_id, nih_counts$gene_id)
  df <- df %>% filter(gene_id %in% common_genes) %>% arrange(gene_id)
  nih <- nih_counts %>% filter(gene_id %in% common_genes) %>% arrange(gene_id)
  mat1 <- as.matrix(df[,-1])
  mat2 <- as.matrix(nih[,-1])
  # Use RMSE instead of mean squared difference
  count_rmse <- sqrt(mean((mat1 - mat2)^2))
  config <- extract_config_count(basename(file))
  tibble(
    count_sd = count_rmse,
    aligner = ifelse(str_detect(file, "star"), "STAR", "HISAT2"),
    min_phred = config$min_phred,
    min_length = config$min_length,
    trim_poly_g = config$trim_poly_g,
    trim_poly_x = config$trim_poly_x
  )
}

summarize_de <- function(file, nih_de) {
  df <- read_csv(file, show_col_types = FALSE)
  colnames(df)[1] <- "gene"
  df <- df %>%
    rename(
      p_value = PValue,
      effect_size = logFC
    )
  common_genes <- intersect(df$gene, nih_de$gene)
  if (length(common_genes) == 0) {
    message("No common genes found in: ", file)
    return(tibble(
      p_value_sd = NA,
      effect_size_sd = NA,
      aligner = ifelse(str_detect(file, "star"), "STAR", "HISAT2"),
      min_phred = NA,
      min_length = NA,
      trim_poly_g = NA,
      trim_poly_x = NA,
      norm_method = NA
    ))
  }
  df <- df %>% filter(gene %in% common_genes) %>% arrange(gene)
  nih <- nih_de %>% filter(gene %in% common_genes) %>% arrange(gene)
  p_value_sd <- mean((df$p_value - nih$p_value)^2)
  effect_size_sd <- mean((df$effect_size - nih$effect_size)^2)
  config <- extract_config_de(basename(file))
  tibble(
    p_value_sd = p_value_sd,
    effect_size_sd = effect_size_sd,
    aligner = ifelse(str_detect(file, "star"), "STAR", "HISAT2"),
    min_phred = config$min_phred,
    min_length = config$min_length,
    trim_poly_g = config$trim_poly_g,
    trim_poly_x = config$trim_poly_x,
    norm_method = config$norm_method
  )
}

# ---- Process all STAR/HISAT2 count matrices ----
star_count_files <- list.files(star_count_dir, pattern = "_matrix.csv$", full.names = TRUE)
hisat2_count_files <- list.files(hisat2_count_dir, pattern = "_matrix.csv$", full.names = TRUE)
count_summaries <- bind_rows(
  lapply(star_count_files, summarize_counts, nih_counts = nih_counts),
  lapply(hisat2_count_files, summarize_counts, nih_counts = nih_counts)
)

# ---- Process all STAR/HISAT2 DE results ----
star_de_files <- list.files(star_de_dir, pattern = "^DE_star_.*\\.csv$", full.names = TRUE)
hisat2_de_files <- list.files(hisat2_de_dir, pattern = "^DE_hisat2_.*\\.csv$", full.names = TRUE)
de_summaries <- bind_rows(
  lapply(star_de_files, summarize_de, nih_de = nih_de),
  lapply(hisat2_de_files, summarize_de, nih_de = nih_de)
)

# ---- Write summary data frames ----
write_csv(count_summaries, "count_sd_df_star_hisat2.csv")
write_csv(de_summaries, "DE_sd_df_star_hisat2.csv")
