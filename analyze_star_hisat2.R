# Libraries
library(dplyr)
library(tidyr)
library(readr)

# Function to process each CSV
process_csv <- function(file_path, nih_count_vector, aligner) {
  # Read the data
  sample_count_df <- read.csv(file_path)
  
  # Extract parameters from the filename
  file_name <- basename(file_path)
  params <- strsplit(file_name, "_")[[1]]
  min_phred <- as.numeric(gsub("Q", "", params[1]))   # Quality threshold (Q22 -> 22)
  min_length <- as.numeric(gsub("L", "", params[2]))  # Length threshold (L35 -> 35)
  trim_poly_g <- as.numeric(params[3] == "G1")        # Trim Poly-G (TRUE -> 1, FALSE -> 0)
  trim_poly_x <- as.numeric(params[4] == "X1")        # Trim Poly-X (TRUE -> 1, FALSE -> 0)
  
  # Reshape the sample counts to long format for easier comparison
  sample_count_long <- sample_count_df %>%
    pivot_longer(cols = -gene_id, names_to = "sample", values_to = "count")
  
  # Join with the NIH counts
  joined_counts <- full_join(
    sample_count_long, nih_count_vector,
    by = c("gene_id" = "gene_id", "sample" = "sample")
  ) %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0)))
  
  # Calculate count standard deviation
  count_var <- sum((joined_counts$count.x - joined_counts$count.y)^2)
  count_sd <- sqrt(count_var / nrow(joined_counts))
  
  # Calculate gene overlap
  overlap_count <- inner_join(
    sample_count_long, nih_count_vector,
    by = c("gene_id" = "gene_id", "sample" = "sample")
  ) %>% nrow()
  overlap_count_percent <- 100 * overlap_count / nrow(joined_counts)
  
  # Collect results
  result <- tibble(
    count_sd = count_sd,
    aligner = aligner,
    min_phred = min_phred,
    min_length = min_length,
    trim_poly_g = trim_poly_g,
    trim_poly_x = trim_poly_x,
    runtime_sec = NA,  # Set to NA for now, as you mentioned
    gene_overlap_percent = overlap_count_percent
  )
  
  return(result)
}

# Main function to iterate over files in the two folders
process_all_files <- function(star_folder, hisat2_folder, nih_count_vector) {
  # List files in each folder
  star_files <- list.files(star_folder, full.names = TRUE)
  hisat2_files <- list.files(hisat2_folder, full.names = TRUE)
  
  # Process each file in STAR folder
  star_results <- lapply(star_files, process_csv, nih_count_vector = nih_count_vector, aligner = "STAR")
  
  # Process each file in HISAT2 folder
  hisat2_results <- lapply(hisat2_files, process_csv, nih_count_vector = nih_count_vector, aligner = "HISAT2")
  
  # Combine all results
  all_results <- bind_rows(star_results, hisat2_results)
  
  # Return the results as a data frame
  return(all_results)
}

# Example usage
nih_count_vector <- read.csv('/lustre/work/client/projects/dheitjan/hschuler_dissertatio/main_storage/STAT-6358-Project/data/gen_samples/nih_counts.csv')
  # Your reference data
colnames(nih_count_vector)[1] = 'gene_id'


star_folder <- "./star_combined_counts"           # Folder with STAR results
hisat2_folder <- "./hisat2_combined_counts"       # Folder with HISAT2 results

# Process all files and get the result
final_results <- process_all_files(star_folder, hisat2_folder, nih_count_vector)

# Output the results to a CSV file
write.csv(final_results, "./STAR_HISAT2_combined_results.csv", row.names = FALSE)
