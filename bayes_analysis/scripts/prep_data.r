# Preps data for Bayesian analysis pipeline.

# Meta Data
library(tidymodels)
tidymodels_prefer()

setwd("/repo/bayes_analysis/")
print(getwd())

sample_names <- c(
    "gene",
    "SRR31476642",
    "SRR31476643",
    "SRR31476644",
    "SRR31476645",
    "SRR31476646",
    "SRR31476647",
    "SRR31476648",
    "SRR31476649",
    "SRR31476650"
)

treatments <- c(
    "DMSO",
    "DMSO",
    "DMSO",
    "DMSO",
    "EPZ015666",
    "EPZ015666",
    "EPZ015666",
    "DMSO",
    "DMSO"
)

factors <- c("aligner", "trim_poly_g", "trim_poly_x", "norm_method")

# Prep count data.
count_sd_df_salmon_kallisto <- read.csv("../data/gen_samples/count_sd_df.csv")
count_sd_df_STAR_HISAT2 <- read.csv("../STAR_HISAT2_combined_results.csv")

count_sd_df <- rbind(count_sd_df_salmon_kallisto, count_sd_df_STAR_HISAT2)
count_sd_df <- count_sd_df |>
    mutate(across(any_of(factors), ~ as.factor(.)))

count_recipe <- count_sd_df |>
    select(-c(runtime_sec, gene_overlap_percent)) |>
    recipe(count_sd ~ .) |>
    step_dummy(all_factor_predictors()) |>
    step_interact(terms = ~ (all_predictors())^2) |>
    # Removes interactions between different levels of the same dummy.
    step_zv(all_predictors()) |>
    step_intercept() |>
    prep()

count_sd_df_clean <- count_recipe |> bake(new_data = NULL)

cat("Writing count sd data to csv... \n")
write.csv(
    count_sd_df_clean, "./data/count_sd_df_clean.csv",
    row.names = FALSE
)

# Prep differential expression data.

DE_sd_df_salmon_kallisto <- read.csv("../data/gen_samples/DE_sd_df.csv")
DE_sd_df <- DE_sd_df_salmon_kallisto
DE_sd_df <- DE_sd_df |>
    mutate(across(any_of(factors), ~ as.factor(.)))

p_value_recipe <- DE_sd_df |>
    select(-c(runtime_sec, gene_overlap_percent, effect_size_sd)) |>
    recipe(p_value_sd ~ .) |>
    step_dummy(all_factor_predictors()) |>
    step_interact(terms = ~ (all_predictors())^2) |>
    # Removes interactions between different levels of the same dummy.
    step_zv(all_predictors()) |>
    step_intercept() |>
    prep()

effect_size_recipe <- DE_sd_df |>
    filter(norm_method != "ALDEx2") |> # Not the same scale as the others.
    droplevels() |>
    select(-c(runtime_sec, gene_overlap_percent, p_value_sd)) |>
    recipe(effect_size_sd ~ .) |>
    step_dummy(all_factor_predictors()) |>
    step_interact(terms = ~ (all_predictors())^2) |>
    # Removes interactions between different levels of the same dummy.
    step_zv(all_predictors()) |>
    step_intercept() |>
    prep()

p_value_sd_df_clean <- p_value_recipe |> bake(new_data = NULL)
effect_size_sd_df_clean <- effect_size_recipe |> bake(new_data = NULL)

cat("Writing p-value data to csv... \n")
write.csv(
    p_value_sd_df_clean, "./data/p_value_sd_df_clean.csv",
    row.names = FALSE
)

cat("Writing effect size data to csv... \n")
write.csv(
    effect_size_sd_df_clean, "./data/effect_size_sd_df_clean.csv",
    row.names = FALSE
)
