# Load libraries
library(tidyverse)

# 1. Import data ##################################################################################################################################

# Access command-line arguments in R
args <- base::commandArgs(trailingOnly = TRUE)

# Assign command-line arguments to R variables
gdepth_file <- args[1]
mismap_file <- args[2]
factor_max_gdepth <- args[3] %>% as.numeric()
median_gdepth_by_individual_0_mismap_file <- args[4]
plot_gdepth_by_individual_file <- args[5]
plot_gdepth_by_region_file <- args[6]
plot_gdepth_by_individual_by_region_file <- args[7]

# Import tab-separated information about genotype depth per site per sample across target sequences presumed to be in single-copy
gdepth <- utils::read.delim(file = gdepth_file, check.names = FALSE)
mismap <- utils::read.delim(file = mismap_file, check.names = FALSE)


# 2. Parse data ###################################################################################################################################

#     2.1 Merge gdepth and mismap =================================================================================================================

# Use pivot longer as very convenient format for merging gdepth and mismap values
# Remove unnecessary columns
gdepth <- gdepth %>%
  dplyr::select(-c(chrom, start, end)) %>%
  tidyr::pivot_longer(cols = 3:(ncol(gdepth)-4),
                      names_to = "sample_id",
                      values_to = "gdepth") %>%
  tidyr::separate_wider_delim(cols = sample_id,
                              delim = "-",
                              names = c("individual_id", "tissue_type")) %>%
  dplyr::select(-tissue_type)
mismap <- mismap %>%
  dplyr::select(-c(chrom, start, end, mismapping_control_targetseq_mean)) %>%
  tidyr::pivot_longer(cols = 3:(ncol(mismap)-5),
                      names_to = "individual_id",
                      values_to = "mismap")

# Merge gdepth and mismap values
gdepth_mismap <- gdepth %>%
  dplyr::full_join(mismap)

#     2.2 Get statistics and plots ================================================================================================================

# Get median by individual
median_gdepth_by_individual <- gdepth_mismap %>%
  dplyr::group_by(individual_id) %>%
  dplyr::summarise(median_gdepth = median(gdepth, na.rm = TRUE),
                   max_gdepth = factor_max_gdepth*median_gdepth)
median_gdepth_by_individual_0_mismap <- gdepth_mismap %>%
  dplyr::filter(mismap == 0) %>%
  dplyr::group_by(individual_id) %>%
  dplyr::summarise(median_gdepth = median(gdepth, na.rm = TRUE),
                   max_gdepth = factor_max_gdepth*median_gdepth)

# Get median by region
median_gdepth_by_region <- gdepth_mismap %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(median_gdepth = median(gdepth, na.rm = TRUE),
                   max_gdepth = factor_max_gdepth*median_gdepth)
median_gdepth_by_region_0_mismap <- gdepth_mismap %>%
  dplyr::filter(mismap == 0) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(median_gdepth = median(gdepth, na.rm = TRUE),
                   max_gdepth = factor_max_gdepth*median_gdepth)

# Get median by individual by region
median_gdepth_by_individual_by_region <- gdepth_mismap %>%
  dplyr::group_by(individual_id, region) %>%
  dplyr::summarise(median_gdepth = median(gdepth, na.rm = TRUE),
                   max_gdepth = factor_max_gdepth*median_gdepth)
median_gdepth_by_individual_by_region_0_mismap <- gdepth_mismap %>%
  dplyr::filter(mismap == 0) %>%
  dplyr::group_by(individual_id, region) %>%
  dplyr::summarise(median_gdepth = median(gdepth, na.rm = TRUE),
                   max_gdepth = factor_max_gdepth*median_gdepth)


# Get histograms

# Plot by individual
plot_gdepth_by_individual <- gdepth_mismap %>%
  ggplot2::ggplot(ggplot2::aes(x = gdepth, fill = !(mismap==0))) +
  ggplot2::geom_histogram(position = "stack", color = "black") +
  ggplot2::scale_fill_manual(values = c("white", "red"), name = "mismap") +
  ggplot2::facet_wrap(~ individual_id, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "gdepth distribution by individual", x = "gdepth") +
  ggplot2::geom_vline(data = median_gdepth_by_individual,
                      ggplot2::aes(xintercept = median_gdepth))

# Plot by region
plot_gdepth_by_region <- gdepth_mismap %>%
  ggplot2::ggplot(ggplot2::aes(x = gdepth, fill = !(mismap==0))) +
  ggplot2::geom_histogram(position = "stack", color = "black") +
  ggplot2::scale_fill_manual(values = c("white", "red"), name = "mismap") +
  ggplot2::facet_wrap(~ region, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "gdepth distribution by regions", x = "gdepth") +
  ggplot2::geom_vline(data = median_gdepth_by_region,
                      ggplot2::aes(xintercept = median_gdepth))

# Plot by individual by region
plot_gdepth_by_individual_by_region <- gdepth_mismap %>%
  ggplot2::ggplot(ggplot2::aes(x = gdepth, fill = !(mismap==0))) +
  ggplot2::geom_histogram(position = "stack", color = "black") +
  ggplot2::scale_fill_manual(values = c("white", "red"), name = "mismap") +
  ggplot2::facet_wrap(individual_id ~ region, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "gdepth distribution by individual by regions", x = "gdepth") +
  ggplot2::geom_vline(data = median_gdepth_by_individual_by_region,
                      ggplot2::aes(xintercept = median_gdepth))


# 3. Export data ##################################################################################################################################

# Save the results to a file

# Median
readr::write_tsv(median_gdepth_by_individual_0_mismap, file = median_gdepth_by_individual_0_mismap_file)

# Plots
ggplot2::ggsave(plot_gdepth_by_individual_file, plot = plot_gdepth_by_individual, width = 12, height = 8)
ggplot2::ggsave(plot_gdepth_by_region_file, plot = plot_gdepth_by_region, width = 12, height = 8)
ggplot2::ggsave(plot_gdepth_by_individual_by_region_file, plot = plot_gdepth_by_individual_by_region, width = 12, height = 8)
