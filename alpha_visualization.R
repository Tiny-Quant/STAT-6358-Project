library(tidyverse)

# Define the output directory
output_dir <- "/users/hschuler/art_hunter/STAT-6358-Project/alpha_output"

# Get a list of all subdirectories
subdirs <- list.dirs(output_dir, recursive = FALSE)

# Initialize empty dataframes
quality_results <- tibble()
alignment_results <- tibble()

# Loop through each subdirectory
for (subdir in subdirs) {
  # Extract metadata from folder name (threshold, length, aligner)
  subdir_name <- basename(subdir)
  match <- str_match(subdir_name, "thr(\\d+)_len(\\d+)_(\\w+)")
  
  if (!is.na(match[1])) {
    threshold <- as.numeric(match[2])
    min_length <- as.numeric(match[3])
    aligner <- match[4]
  } else {
    next  # Skip directories that don't match the expected naming pattern
  }
  
  # Process quality_log.csv
  quality_file <- file.path(subdir, "quality_log.csv")
  if (file.exists(quality_file)) {
    q_data <- read_csv(quality_file, col_names = c("filename", "threshold", "min_length", "input_reads", "output_reads", "processing_time"), show_col_types = FALSE)
    q_data <- q_data %>% mutate(aligner = aligner)
    quality_results <- bind_rows(quality_results, q_data)
  }
  
  # Process alignment_log.csv
  alignment_file <- file.path(subdir, "alignment_log.csv")
  if (file.exists(alignment_file) && file.info(alignment_file)$size > 0) {  # Skip empty files
    a_data <- read_csv(alignment_file, col_names = c("filename", "aligner", "aligned_reads"), show_col_types = FALSE)
    a_data <- a_data %>% mutate(aligned_reads = as.numeric(aligned_reads))  # Ensure 'aligned_reads' is numeric
    alignment_results <- bind_rows(alignment_results, a_data)
  }
}

# Merge results on filename and aligner
final_results <- full_join(quality_results, alignment_results, by = c("filename", "aligner"))

# Save combined results as a CSV
write_csv(final_results, file.path(output_dir, "combined_results.csv"))

# Print summary
print(head(final_results))

# Generate quality metric plots
ggplot(final_results, aes(x = factor(threshold), y = output_reads, fill = aligner)) +
  geom_boxplot() +
  facet_wrap(~ min_length) +
  labs(title = "Output Reads Across Thresholds and Aligners", x = "Quality Threshold", y = "Output Reads") +
  theme_minimal() +
  ggsave(filename = file.path(output_dir, "output_reads_plot.png"), width = 10, height = 6)

ggplot(final_results, aes(x = factor(threshold), y = processing_time, fill = aligner)) +
  geom_boxplot() +
  facet_wrap(~ min_length) +
  labs(title = "Processing Time Across Thresholds and Aligners", x = "Quality Threshold", y = "Processing Time (s)") +
  theme_minimal() +
  ggsave(filename = file.path(output_dir, "processing_time_plot.png"), width = 10, height = 6)

# Generate alignment rate plot (only for non-missing values)
if ("aligned_reads" %in% names(final_results)) {
  ggplot(final_results %>% drop_na(aligned_reads), aes(x = factor(threshold), y = aligned_reads, fill = aligner)) +
    geom_boxplot() +
    facet_wrap(~ min_length) +
    labs(title = "Aligned Reads Across Thresholds and Aligners", x = "Quality Threshold", y = "Aligned Reads") +
    theme_minimal() +
    ggsave(filename = file.path(output_dir, "aligned_reads_plot.png"), width = 10, height = 6)
}

cat("Analysis complete. Results saved in", output_dir, "\n")
