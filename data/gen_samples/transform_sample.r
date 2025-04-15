# Libraries
library(tidymodels)
library(edgeR)
library(tximport)

# Create baseline files if they don't already exist. 

# Meta-Data from:  
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1189593&o=acc_s%3Aa

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

if (file.exists("nih_counts.csv")) {
    message("nih_counts.csv already exists.")
} else {
    message("Generating: nih_counts.csv")

    nih_count_matrix <- read.csv("../GSE282674_ovcar4_count_table.csv")
    colnames(nih_count_matrix) <- sample_names

    nih_count_matrix <- nih_count_matrix |>
        mutate(gene = str_remove(gene, "\\..*$"))

    nih_count_vector <- nih_count_matrix |> 
        pivot_longer(cols = -gene, names_to = "sample", values_to = "count")

    write.csv(nih_count_vector, file = "nih_counts.csv", row.names = FALSE)
}

if (file.exists("nih_p_values.csv")) {
    message("nih_p_values.csv already exists.")
} else {
    message("Generating: nih_p_values.csv")

    nih_count_matrix <- read.csv("../GSE282674_ovcar4_count_table.csv")
    colnames(nih_count_matrix) <- sample_names
    nih_count_matrix <- nih_count_matrix |>
        mutate(gene = str_remove(gene, "\\..*$"))

    nih_dgelist <- DGEList(
        counts = nih_count_matrix[, -1], 
        genes = nih_count_matrix$gene, 
        group = as.factor(treatments)
    )
    nih_dgelist <- nih_dgelist[filterByExpr(nih_dgelist), , keep.lib.sizes = FALSE]
    nih_dgelist <- estimateDisp(nih_dgelist)

    design <- model.matrix(~as.factor(treatments))
    nih_fit <- glmFit(nih_dgelist, design)
    nih_LRT <- glmLRT(nih_fit, coef = 2)

    nih_p_values <- topTags(nih_LRT, n = Inf)$table
    nih_p_values <- data.frame(
        gene = nih_p_values$genes, 
        p_value = nih_p_values$PValue
    )

    write.csv(nih_p_values, file = "nih_p_values.csv", row.names = FALSE)
}

# Read in sample count matrix from pipeline. 
args <- commandArgs(trailingOnly = TRUE)

sample_dir   <- args[1]
aligner       <- args[2]
min_phred     <- args[3]
min_length    <- args[4]
trim_poly_g   <- args[5]
trim_poly_x   <- args[6]
runtime_sec   <- as.numeric(args[7])

# Generate count matrix. 
tx2gene <- read.csv("../align_indices/tx2gene.csv")
sample_dirs <- list.dirs(sample_dir, recursive = FALSE, full.name = TRUE)

if (aligner == "salmon") {
   files <- file.path(sample_dirs, "quant.sf") 
   names(files) <- basename(sample_dirs)
   txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                   ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)

    sample_count_df <- txi$counts |> as.data.frame() |>
                       rownames_to_column(var = "X") 
    colnames(sample_count_df) <- sample_names

} else if(aligner == "kallisto") {
    files <- file.path(sample_dirs, "abundance.h5")
    names(files) <- basename(sample_dirs)
    txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, 
                    ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)

    sample_count_df <- txi$counts |> as.data.frame() |>
                       rownames_to_column(var = "X") 
    colnames(sample_count_df) <- sample_names

} else if (aligner == "rsem") {
    files <- list.files(
        sample_dir, pattern = ".genes.results$", full.names = TRUE
    )
    names(files) <- basename(sample_dirs)
    print(sample_dir)
    txi <- tximport(files, type = "rsem", tx2gene = tx2gene, 
                    txIn = FALSE, txOut = FALSE, 
                    ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)

    sample_count_df <- txi$counts |> as.data.frame() |>
                       rownames_to_column(var = "X") 
    colnames(sample_count_df) <- sample_names

} else {
    message("Invalid aligner: ", aligner)
}

# Compute count statistic. 

nih_count_vector <- read.csv("./nih_counts.csv")

sample_count_long <- sample_count_df |> 
    pivot_longer(cols = -gene, names_to = "sample", values_to = "count")

joined_counts <- full_join(
    sample_count_long, nih_count_vector, 
    by = join_by(gene, sample) 
) |> mutate(across(where(is.numeric), ~replace_na(., 0)))

count_var <- sum((joined_counts$count.x - joined_counts$count.y)^2)
count_sd <- sqrt(count_var / nrow(joined_counts))

# Calculate overlap. 
overlap_count <- inner_join(
    sample_count_long, nih_count_vector, 
    by = join_by(gene, sample) 
) |> nrow()
overlap_count_percent <- 100 * overlap_count / nrow(joined_counts)


sample_count_data <- c(
    count_sd, 
    aligner, min_phred, min_length, trim_poly_g, trim_poly_x, 
    runtime_sec, overlap_count_percent 
)

# Compute p-value statistic. 

sample_dgelist <- DGEList(
    counts = sample_count_df[, -1], 
    genes = sample_count_df$gene, 
    group = as.factor(treatments)
)

nih_p_values <- read.csv("./nih_p_values.csv")

sample_dgelist <- sample_dgelist[
    filterByExpr(sample_dgelist), , keep.lib.sizes = FALSE
]
sample_dgelist <- estimateDisp(sample_dgelist)

# Random select a normalization method. 
norm_methods <- c("default", "TMM", "RLE", "upperquartile")
method_choice <- sample(norm_methods, 1)

if (method_choice == "default") {
    sample_dgelist <- calcNormFactors(sample_dgelist)
} else {
    sample_dgelist <- calcNormFactors(sample_dgelist, method = method_choice)
}

design <- model.matrix(~as.factor(treatments))
sample_fit <- glmFit(sample_dgelist, design)
sample_LRT <- glmLRT(sample_fit, coef = 2)

sample_p_values <- topTags(sample_LRT, n = Inf)$table
sample_p_values <- data.frame(
    gene = sample_p_values$genes, 
    p_value = sample_p_values$PValue
)

joined_p_values <- full_join(
    sample_p_values, nih_p_values, by = join_by(gene)
) |> mutate(across(where(is.numeric), ~replace_na(., 1)))

p_value_var <- sum((joined_p_values$p_value.x - joined_p_values$p_value.y)^2)
p_value_sd <- sqrt(p_value_var / nrow(joined_p_values))

# Calculate overlap. 
overlap_p_value <- inner_join(
    sample_p_values, nih_p_values, by = join_by(gene)
) |> nrow()
overlap_p_value_percent <- 100 * overlap_p_value / nrow(joined_p_values)

sample_p_value_data <- c(
    p_value_sd, 
    aligner, min_phred, min_length, trim_poly_g, trim_poly_x, 
    runtime_sec, 
    method_choice, overlap_p_value_percent
)

# Output sample. 
unique_id <- paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample(10000:99999, 1)
)

write.csv(
    sample_count_data, 
    file = paste0(
        "./count_sd_samples/", "count_sd_sample_", unique_id, ".csv"
    ), 
    row.names = FALSE 
)

write.csv(
    sample_p_value_data, 
    file = paste0(
        "./p_value_sd_samples/", "p_value_sd_sample_", unique_id, ".csv"
    ), 
    row.names = FALSE 
)