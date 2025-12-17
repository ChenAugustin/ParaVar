# Load libraries
library(tidyverse)

# 1. Import data ##################################################################################################################################

# Access command-line arguments in R
args <- base::commandArgs(trailingOnly = TRUE)

# Assign command-line arguments to R variables
gdepth_controllib_targetseq_file <- args[1]
gdepth_targetlib_targetseq_file <- args[2]
library_size_file <- args[3]
fraction_targetseq_in_targetlib <- args[4] %>% as.numeric()
output_tsv_file <- args[5]

# Import comma-separated information about library size per sample
library_size <- utils::read.csv(file = library_size_file)

# Import tab-separated information about genotype depth per site per sample across target sequence
gdepth_controllib_targetseq <- utils::read.delim(file = gdepth_controllib_targetseq_file)
gdepth_targetlib_targetseq <- utils::read.delim(file = gdepth_targetlib_targetseq_file)


# 2. Parse data ###################################################################################################################################

#     2.1 Get library size ratio ==================================================================================================================

# Compute library size ratio per individual between targetlib and controllib across control sequence:
#     1. Select relevant columns
#     2. Remove duplicates by keeping one row per sample ID
#     3. Edit the "targetlib" column for it to contains the suffixes of the new "library_size" columns after pivot_wider()
#     4. Remove sample ID as impedes merging based on individual ID
#     5. Have matched controllib and targetlib library sizes on 2 different columns using pivot_wider()
#     6. Which allows to calculate the library size ratio between the matched targetlib and controllib for each individual ID
#     7. Remove library sizes and keep only library size ratios to lighten the data frame
library_size_ratios_controlseq <- library_size %>%
  dplyr::select(targetlib, individual_id, sample_id, library_size) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(targetlib = dplyr::if_else(targetlib == "TRUE", "targetlib", "controllib")) %>%
  dplyr::select(-sample_id) %>%
  tidyr::pivot_wider(names_from = targetlib, values_from = library_size, names_prefix = "library_size_") %>%
  dplyr::mutate(library_size_ratio_controlseq = (library_size_targetlib * (1 - fraction_targetseq_in_targetlib)) / library_size_controllib) %>%
  dplyr::select(-c(library_size_targetlib, library_size_controllib))

#     2.2 Prepare genotype depth info for merging =================================================================================================

# Create function to:
#     1. Pivot genotype depth data to long format to reveal sites with missing genotype depth value after merging
#     2. Extract individual_id from sample_id for merging
#     3. Remove sample_id as impedes merging based on individual_id
pivot_longer_gdepth <- function(df, value_name) {
  df %>%
    tidyr::pivot_longer(cols = -c(CHROM, POS), names_to = "sample_id", values_to = value_name) %>%
    dplyr::mutate(individual_id = stringr::str_split_fixed(sample_id, "\\.", 2)[,1]) %>%
    dplyr::select(-sample_id)
}

# Use function to pivot controllib and targetlib genotype depth data to long format
# Name the gdepth column of each file "gdepth_controllib_targetseq" and "gdepth_targetlib_targetseq" to distinguish them cleanly when merging
gdepth_controllib_targetseq_long <- pivot_longer_gdepth(gdepth_controllib_targetseq, "gdepth_controllib_targetseq")
gdepth_targetlib_targetseq_long <- pivot_longer_gdepth(gdepth_targetlib_targetseq, "gdepth_targetlib_targetseq")

#     2.3 Merge all data and get mismapping of control reads on targetseq =========================================================================

# Merge all data:
#     1. Merge both genotype depth data frames with full_join() to reveal sites with missing genotype depth value
#     2. Re-sort by chromosome name and genomic position as it gets mixed up at this step
#     3. Fill missing genotype depth values with 0
#     4. And join with library size ratios across controlseq using individual ID
# Compute mismapping of control reads on targetseq:
#     1. Calculate mismapping of control reads on targetseq per site per individual
#     2. Replace NaN by 0 and Inf by 1
#     2. Remove genotype depth as irrelevant here
#     3. Get final desired format: chromosome, genomic position, and mismapping of control reads on targetseq per individual per site
#     4. Add mean mismapping of control reads on targetseq per site across all individuals
mismapping_data_control_targetseq <- dplyr::full_join(gdepth_controllib_targetseq_long, gdepth_targetlib_targetseq_long) %>%
  dplyr::arrange(CHROM, POS) %>%
  tidyr::replace_na(list(gdepth_controllib_targetseq = 0, gdepth_targetlib_targetseq = 0)) %>%
  dplyr::left_join(library_size_ratios_controlseq, by = "individual_id") %>%
  dplyr::mutate(mismapping_control_targetseq = (gdepth_controllib_targetseq/gdepth_targetlib_targetseq) * library_size_ratio_controlseq) %>%
  tidyr::replace_na(list(mismapping_control_targetseq = 0)) %>%
  dplyr::mutate(mismapping_control_targetseq = ifelse(is.infinite(mismapping_control_targetseq), 1, mismapping_control_targetseq)) %>%
  dplyr::select(CHROM, POS, individual_id, mismapping_control_targetseq) %>%
  tidyr::pivot_wider(names_from = individual_id, values_from = mismapping_control_targetseq) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(mismapping_control_targetseq_mean = base::mean(dplyr::c_across(-c(CHROM, POS)), na.rm = TRUE))


# 3. Export data ##################################################################################################################################

# Save the result to a file
readr::write_tsv(mismapping_data_control_targetseq, file = output_tsv_file)
