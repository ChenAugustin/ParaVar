# Load libraries
library(tidyverse)
library(RIdeogram)

# 1. Import data ##################################################################################################################################

# Access command-line arguments in R
args <- base::commandArgs(trailingOnly = TRUE)

# Assign command-line arguments to R variables
targetseq_file <- args[1]
pixy_bed_file <- args[2]
pixy_output_file <- args[3]
plot_scorability_file <- args[4]
plot_pi_file <- args[5]

# Import tab-separated information about windows and statistics per window
targetseq <- utils::read.delim(file = targetseq_file, col.names = c("Chr", "Start", "End"))
pixy_bed <- utils::read.delim(file = pixy_bed_file, col.names = c("Chr", "Start", "End"))
pixy_output <- utils::read.delim(file = pixy_output_file)


# 2. Parse data ###################################################################################################################################

#     2.1 Extract statistics ======================================================================================================================

# Extract statistics from pixy's output
pi <- pixy_output %>%
  dplyr::select(Chr=chromosome, Start=window_pos_1, End=window_pos_2, Value=avg_pi)
scorability <- pixy_output %>%
  dplyr::select(Chr=chromosome, Start=window_pos_1, End=window_pos_2, Value=no_sites)

# Add contigs missing from pixy output because 0 scorable site
pi <- pi %>% dplyr::full_join(pixy_bed)
scorability <- scorability %>% dplyr::full_join(pixy_bed)

# Contigs with 0 scorability end up with NA scorability for all Windows, an artefact from full_join()
# Replace NA by 0 as it should be
scorability <- scorability %>%
  dplyr::mutate(Value = dplyr::if_else(is.na(Value), 0, Value))

#     2.2 Plot ====================================================================================================================================

# Plot scorability along all 692 contigs
#     1. Convert objects to data.frame (instead of tibbles for example)
#     2. Order the contigs from the longest to the smallest. The contigs will appear in this order
#     3. Select color gradient
#     4. Define output name, size, and save as .svg
RIdeogram::ideogram(karyotype = as.data.frame(targetseq %>% dplyr::arrange(desc(End))),
                    overlaid = scorability %>% as.data.frame(),
                    colorset1 = c("#4575b4", "#ffffbf", "#cc6677"),
                    output = plot_scorability_file,
                    width = 1700)

# Plot pi along all target sequences
#     1. Same process as for scorability
RIdeogram::ideogram(karyotype = as.data.frame(targetseq %>% dplyr::arrange(desc(End))),
                    overlaid = pi %>% as.data.frame(),
                    colorset1 = c("#FFFFFF", "#FD00FD"),
                    output = plot_pi_file,
                    width = 1700)
