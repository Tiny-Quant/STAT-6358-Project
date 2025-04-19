# Libraries
library(tidymodels)
library(edgeR)
library(tximport)
library(ALDEx2)

#############
# Meta-Data #
#############
# ref: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1189593&o=acc_s%3Aa

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

#########################
# Create Baseline Files #
#########################

if (file.exists("nih_counts.csv")) {
    message("nih_counts.csv already exists.")
} else {
    message("Generating: nih_counts.csv")

    nih_count_matrix <- read.csv("../GSE282674_ovcar4_count_table.csv")
    colnames(nih_count_matrix) <- sample_names

    nih_count_matrix <- nih_count_matrix |>
        mutate(gene = stringr::str_remove(gene, "\\..*$"))

    nih_count_vector <- nih_count_matrix |>
        pivot_longer(cols = -gene, names_to = "sample", values_to = "count")

    write.csv(nih_count_vector, file = "nih_counts.csv", row.names = FALSE)
}

if (file.exists("nih_DE.csv")) {
    message("nih_DE.csv already exists.")
} else {
    message("Generating: nih_DE.csv")

    nih_count_matrix <- read.csv("../GSE282674_ovcar4_count_table.csv")
    colnames(nih_count_matrix) <- sample_names
    nih_count_matrix <- nih_count_matrix |>
        mutate(gene = stringr::str_remove(gene, "\\..*$"))

    nih_dgelist <- DGEList(
        counts = nih_count_matrix[, -1],
        genes = nih_count_matrix$gene,
        group = as.factor(treatments)
    )

    nih_dgelist <- nih_dgelist[
        filterByExpr(nih_dgelist), ,
        keep.lib.sizes = FALSE
    ]

    nih_dgelist <- estimateDisp(nih_dgelist)

    nih_dgelist <- calcNormFactors(nih_dgelist, method = "none")

    design <- model.matrix(~ as.factor(treatments))
    nih_fit <- glmFit(nih_dgelist, design)
    nih_LRT <- glmLRT(nih_fit, coef = 2)

    nih_DE <- topTags(nih_LRT, n = Inf)$table
    nih_DE <- data.frame(
        gene = nih_DE$genes,
        p_value = nih_DE$PValue,
        effect_size = nih_DE$logFC
    )

    write.csv(nih_DE, file = "nih_DE.csv", row.names = FALSE)
}

################################################################################

start_time <- proc.time()[3]

# Read in sample count matrix from pipeline.
args <- commandArgs(trailingOnly = TRUE)

sample_dir <- args[1]
print(sample_dir)
aligner <- args[2]
min_phred <- args[3]
min_length <- args[4]
trim_poly_g <- args[5]
trim_poly_x <- args[6]
runtime_sec <- as.numeric(args[7])

#########################
# Generate Count Matrix #
#########################

tx2gene <- read.csv("../align_indices/tx2gene.csv")
sample_dirs <- list.dirs(sample_dir, recursive = FALSE, full.name = TRUE)

if (aligner == "salmon") {
    valid_dirs <- sample_dirs[
        file.exists(file.path(sample_dirs, "quant.sf"))
    ]
    files <- file.path(valid_dirs, "quant.sf")
    names(files) <- basename(valid_dirs)
    print(files)
    txi <- tximport(files,
        type = "salmon", tx2gene = tx2gene,
        ignoreTxVersion = TRUE, ignoreAfterBar = TRUE
    )

    sample_count_df <- txi$counts |>
        as.data.frame() |>
        rownames_to_column(var = "X")
    colnames(sample_count_df) <- sample_names
} else if (aligner == "kallisto") {
    valid_dirs <- sample_dirs[
        file.exists(file.path(sample_dirs, "abundance.h5"))
    ]
    files <- file.path(valid_dirs, "abundance.h5")
    names(files) <- basename(valid_dirs)
    print(files)
    txi <- tximport(files,
        type = "kallisto", tx2gene = tx2gene,
        ignoreTxVersion = TRUE, ignoreAfterBar = TRUE
    )

    sample_count_df <- txi$counts |>
        as.data.frame() |>
        rownames_to_column(var = "X")
    colnames(sample_count_df) <- sample_names
} else if (aligner == "rsem") {
    # files <- list.files(
    #     sample_dir,
    #     pattern = ".genes.results$", full.names = TRUE
    # )
    # names(files) <- basename(sample_dirs)
    valid_dirs <- sample_dirs[
        file.exists(file.path(sample_dirs, "genes.results"))
    ]
    files <- file.path(valid_dirs, "genes.results")
    names(files) <- basename(valid_dirs)
    print(files)
    txi <- tximport(files,
        type = "rsem", tx2gene = tx2gene,
        txIn = FALSE, txOut = FALSE,
        ignoreTxVersion = TRUE, ignoreAfterBar = TRUE
    )

    sample_count_df <- txi$counts |>
        as.data.frame() |>
        rownames_to_column(var = "X")
    colnames(sample_count_df) <- sample_names
} else {
    message("Invalid aligner: ", aligner)
}

############################
# Compute Count Statistics #
############################

nih_count_vector <- read.csv("./nih_counts.csv")

sample_count_long <- sample_count_df |>
    pivot_longer(cols = -gene, names_to = "sample", values_to = "count")

joined_counts <- full_join(
    sample_count_long, nih_count_vector,
    by = join_by(gene, sample)
) |> mutate(across(where(is.numeric), ~ replace_na(., 0)))

count_var <- sum((joined_counts$count.x - joined_counts$count.y)^2)
count_sd <- sqrt(count_var / nrow(joined_counts))

# Calculate overlap.
overlap_count <- inner_join(
    sample_count_long, nih_count_vector,
    by = join_by(gene, sample)
) |> nrow()
overlap_count_percent <- 100 * overlap_count / nrow(joined_counts)

count_end_time <- proc.time()[3]

sample_count_data <- c(
    count_sd,
    aligner, min_phred, min_length, trim_poly_g, trim_poly_x,
    runtime_sec + count_end_time - start_time,
    overlap_count_percent
)

##############################################
# Compute Differential Expression Statistics #
##############################################

nih_DE <- read.csv("./nih_DE.csv")

# Random select a normalization method.
norm_methods <- c("none", "TMM", "TMMwsp", "RLE", "upperquartile", "ALDEx2")
method_choice <- sample(norm_methods, 1)

if (method_choice == "ALDEx2") {
    aldex_count_matrix <- (
        sample_count_df |>
            column_to_rownames(var = "gene") |>
            mutate(
                across(
                    where(is.numeric), ~ round(as.integer(.))
                )
            )
    )

    aldex_results <- aldex(
        reads = aldex_count_matrix,
        conditions = treatments |> as.factor() |> as.numeric()
    )

    sample_DE <- data.frame(
        gene = rownames(aldex_results),
        p_value = aldex_results$we.eBH,
        effect_size = aldex_results$effect
    )
} else {
    sample_dgelist <- DGEList(
        counts = sample_count_df[, -1],
        genes = sample_count_df$gene,
        group = as.factor(treatments)
    )

    sample_dgelist <- sample_dgelist[
        filterByExpr(sample_dgelist), ,
        keep.lib.sizes = FALSE
    ]

    sample_dgelist <- estimateDisp(sample_dgelist)

    sample_dgelist <- calcNormFactors(sample_dgelist, method = method_choice)

    design <- model.matrix(~ as.factor(treatments))
    sample_fit <- glmFit(sample_dgelist, design)
    sample_LRT <- glmLRT(sample_fit, coef = 2)

    sample_DE <- topTags(sample_LRT, n = Inf)$table
    sample_DE <- data.frame(
        gene = sample_DE$genes,
        p_value = sample_DE$PValue,
        effect_size = sample_DE$logFC
    )
}

joined_DE <- full_join(
    sample_DE, nih_DE,
    by = join_by(gene)
) |> mutate(across(where(is.numeric), ~ replace_na(., 1)))

p_value_var <- sum((joined_DE$p_value.x - joined_DE$p_value.y)^2)
p_value_sd <- sqrt(p_value_var / nrow(joined_DE))

effect_size_var <- sum((joined_DE$effect_size.x - joined_DE$effect_size.y)^2)
effect_size_sd <- sqrt(effect_size_var / nrow(joined_DE))

# Calculate overlap.
overlap_DE <- inner_join(
    sample_DE, nih_DE,
    by = join_by(gene)
) |> nrow()
overlap_DE_percent <- 100 * overlap_DE / nrow(joined_DE)

DE_end_time <- proc.time()[3]

sample_DE_data <- c(
    p_value_sd,
    effect_size_sd,
    aligner, min_phred, min_length, trim_poly_g, trim_poly_x,
    runtime_sec + DE_end_time - count_end_time,
    method_choice, overlap_DE_percent
)

#################
# Output sample #
#################

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
    sample_DE_data,
    file = paste0(
        "./DE_sd_samples/", "DE_sd_sample_", unique_id, ".csv"
    ),
    row.names = FALSE
)
