#############
# Meta Data #
#############
library(tidyverse)

setwd("/repo/data/gen_samples")

factors <- c("aligner", "trim_poly_g", "trim_poly_x", "norm_method")

#####################
# Format Count Data #
#####################
count_sample_dir <- "./count_sd_samples"
count_files <- list.files(
    count_sample_dir,
    pattern = "\\.csv", full.names = TRUE
)
count_list <- lapply(count_files, read.csv)

count_sd_df <- do.call(cbind, count_list) |>
    t() |>
    as.data.frame()
rownames(count_sd_df) <- NULL
colnames(count_sd_df) <- c(
    "count_sd", "aligner", "min_phred", "min_length",
    "trim_poly_g", "trim_poly_x", "runtime_sec",
    "gene_overlap_percent"
)
count_sd_df <- count_sd_df |>
    mutate(across(any_of(factors), ~ as.factor(.))) |>
    mutate(across(-any_of(factors), ~ as.numeric(.)))

write.csv(count_sd_df, "./count_sd_df.csv", row.names = FALSE)

#######################################
# Format Differential Expression Data #
#######################################
DE_sample_dir <- "./DE_sd_samples"
DE_files <- list.files(
    DE_sample_dir,
    pattern = "\\.csv", full.names = TRUE
)
DE_list <- lapply(DE_files, read.csv)

DE_sd_df <- do.call(cbind, DE_list) |>
    t() |>
    as.data.frame()
rownames(DE_sd_df) <- NULL
colnames(DE_sd_df) <- c(
    "p_value_sd", "effect_size_sd", "aligner", "min_phred", "min_length",
    "trim_poly_g", "trim_poly_x", "runtime_sec", "norm_method",
    "gene_overlap_percent"
)
DE_sd_df <- DE_sd_df |>
    mutate(across(any_of(factors), ~ as.factor(.))) |>
    mutate(across(-any_of(factors), ~ as.numeric(.)))

write.csv(DE_sd_df, "./DE_sd_df.csv", row.names = FALSE)
