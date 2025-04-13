library(ggplot2)
library(dplyr)
library(tidyr)

# Define directories
output_dir <- "alpha_output"
thresh_levels <- c(25, 30, 39)

# Function to parse STAR log.final.out
parse_star_log <- function(log_file) {
  log_data <- readLines(log_file)
  extract_value <- function(pattern) {
    line <- log_data[grep(pattern, log_data)]
    if (length(line) > 0) {
      return(as.numeric(strsplit(line, "\t")[[1]][2]))
    }
    return(NA)
  }
  
  total_reads <- extract_value("Number of input reads")
  unique_reads <- extract_value("Uniquely mapped reads %")
  multi_mapped <- extract_value("% of reads mapped to multiple loci")
  return(data.frame(Total=total_reads, Unique=unique_reads, Multi=multi_mapped))
}

# Function to parse RSubread summary
parse_rsubread_log <- function(summary_file) {
  df <- read.delim(summary_file, header=TRUE)
  total <- sum(df$Assigned + df$Unassigned_NoFeatures + df$Unassigned_Ambiguity + 
                 df$Unassigned_MultiMapping + df$Unassigned_Unmapped)
  assigned <- sum(df$Assigned) / total * 100
  return(data.frame(Total=total, Assigned=assigned))
}

# Collect data
alignment_data <- data.frame()

for (thresh in thresh_levels) {
  star_log <- file.path(output_dir, paste0("thr", thresh, "_STAR"), "Log.final.out")
  rsubread_summary <- file.path(output_dir, paste0("thr", thresh, "_Rsubread"), "summary.txt")
  
  if (file.exists(star_log)) {
    star_results <- parse_star_log(star_log)
    star_results$Threshold <- thresh
    star_results$Aligner <- "STAR"
    alignment_data <- rbind(alignment_data, star_results)
  }
  
  if (file.exists(rsubread_summary)) {
    rsubread_results <- parse_rsubread_log(rsubread_summary)
    rsubread_results$Threshold <- thresh
    rsubread_results$Aligner <- "RSubread"
    alignment_data <- rbind(alignment_data, rsubread_results)
  }
}

# Reshape for plotting
alignment_long <- alignment_data %>% 
  pivot_longer(cols = c(Unique, Multi, Assigned), names_to = "Metric", values_to = "Percentage")

# Plot alignment rates
p <- ggplot(alignment_long, aes(x = factor(Threshold), y = Percentage, fill = Aligner)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Alignment Rates by Quality Threshold", x = "Quality Threshold", y = "Percentage (%)")

print(p)
