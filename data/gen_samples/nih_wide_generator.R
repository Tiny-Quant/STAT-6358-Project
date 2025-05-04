library(readr)
library(tidyr)
library(dplyr)

nih_counts_long <- read_csv("/lustre/work/client/projects/dheitjan/hschuler_dissertatio/main_storage/STAT-6358-Project/data/gen_samples/nih_counts.csv")
nih_counts_wide <- nih_counts_long %>%
  pivot_wider(names_from = sample, values_from = count) %>%
  rename(gene_id = gene)

write_csv(nih_counts_wide, "/lustre/work/client/projects/dheitjan/hschuler_dissertatio/main_storage/STAT-6358-Project/data/gen_samples/nih_counts_wide.csv")