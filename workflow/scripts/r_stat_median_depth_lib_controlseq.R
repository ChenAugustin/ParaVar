# Load libraries
library(tidyverse)

# 1. Import data ##################################################################################################################################

# Access command-line arguments in R
args <- base::commandArgs(trailingOnly = TRUE)

# Assign command-line arguments to R variables
shuf_samtools_depth_lib_controlseq_file = args[1]
output_tsv <-  args[2]

# Import tab-separated information about depth per site per sample across a random sample of the control sequence
# Do not edit the column names (e.g., convert "-" to ".") by setting check.names = FALSE
shuf_samtools_depth_lib_controlseq <- utils::read.delim(file = shuf_samtools_depth_lib_controlseq_file, check.names = FALSE)


# 2. Parse data ###################################################################################################################################

# Calculate median depth per sample
shuf_samtools_median_depth_lib_controlseq_per_sample <- shuf_samtools_depth_lib_controlseq %>%
  select(-CHROM, -POS) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "sample_id", values_to = "median_depth_lib_controlseq")


# 3. Export data ##################################################################################################################################

# Save the result to a file
readr::write_tsv(shuf_samtools_median_depth_lib_controlseq_per_sample, file = output_tsv)
