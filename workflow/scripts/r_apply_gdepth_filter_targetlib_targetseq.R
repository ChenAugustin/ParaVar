# Load libraries
library(tidyverse)

# 1. Import data ##################################################################################################################################

# Access command-line arguments in R
args <- base::commandArgs(trailingOnly = TRUE)

# Assign command-line arguments to R variables
gdepth_file <- args[1]
median_gdepth_by_individual_0_mismap_file <- args[2]
min_gdepth <- args[3]
output_gdepth_file <- args[4]

# Import tab-separated information about genotype depth and median gdepth by individual across target single-copy regions
gdepth <- utils::read.delim(file = gdepth_file, check.names = FALSE)
median_gdepth_by_individual_0_mismap <- utils::read.delim(file = median_gdepth_by_individual_0_mismap_file, check.names = FALSE)


# 2. Parse data ###################################################################################################################################

#     2.1 Reformat data ===========================================================================================================================

# The gdepth data

# Convert sample_id to individual_id
base::names(gdepth)[3:base::ncol(gdepth)] <- base::sub("-.*", "", base::names(gdepth)[3:base::ncol(gdepth)])
# Convert into data.table for memory-efficient processing
dt <- data.table::as.data.table(gdepth)
# Select columns with gdepth information
cols <- base::names(dt)[3:base::ncol(gdepth)]
# create a logical copy of just those columns (so original dt remains intact)
dt_logic <- dt[, ..cols]              # shallow copy of references to columns

# The upper gdepth thresholds

# Use wider format
median_gdepth_by_individual_0_mismap <- median_gdepth_by_individual_0_mismap %>%
  dplyr::select(-median_gdepth) %>%
  tidyr::pivot_wider(values_from = max_gdepth,
                     names_from = individual_id)
# Then convert thresholds data frame to vector format
# While rearranging the values to get the same order as in the gdepth data
thr <- base::as.numeric(median_gdepth_by_individual_0_mismap[1, cols, drop = TRUE])

#     2.2 Apply gdepth filter =====================================================================================================================

# Loop over the indices of the selected columns
for (i in base::seq_along(cols)) {
  
  # In the logical data.table (dt_logic), replace column i with a logical vector
  # created by testing the corresponding column in the original data.table (dt).
  # Each element is TRUE if:
  #   - the value in dt[[cols[i]]] is >= min_gdepth
  #   - AND the value is <= the column-specific threshold thr[i]
  #
  # Arguments to data.table::set():
  #   j = i        → column index (which column to replace)
  #   value = ...  → new vector (must be the same length as number of rows)
  #
  # set() modifies dt_logic in place (no copy), so it is memory-efficient.
  data.table::set(
    dt_logic,              # data.table we are modifying
    j = i,                 # which column to update (by numeric index here)
    value = dt[[cols[i]]] >= min_gdepth & dt[[cols[i]]] <= thr[i]  # logical vector
  )
}

# For each site, count the number of individuals that passed the depth filter 
n_depth_pass <- dt_logic[, base::rowSums(.SD), .SDcols = cols]

# Add this information to original gdepth data
gdepth <- base::cbind(gdepth, n_depth_pass)


# 3. Export data ##################################################################################################################################

readr::write_tsv(gdepth, file = output_gdepth_file)
