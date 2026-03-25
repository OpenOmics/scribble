# Load the necessary libraries
library(hdf5r)
library(tidyverse) # Loads tibble, readr, dplyr, etc.

# 1. Define file paths
cellbender_file <- 'Sample_corrected.h5'
output_txt_file <- 'Sample_probability.txt'

# 2. Open the HDF5 file
# R's error handling will automatically stop if the file isn't found
h5_file <- H5File$new(cellbender_file, mode = 'r')

# --- Extract the data ---

  # 1. Read all original barcodes
full_barcodes <- h5_file[['matrix/barcodes']]$read()

# 2. Read the indices that link cell_probability to full_barcodes
# HDF5 indices are 0-based, R indices are 1-based, so add 1
  barcode_indices_for_latents <- h5_file[['droplet_latents/barcode_indices_for_latents']]$read() + 1

  # 3. Use these indices to get the subset of barcodes corresponding to cell_probability
  # These will be the barcodes for which CellBender computed latent variables
  subset_barcodes <- full_barcodes[barcode_indices_for_latents]

  # 4. Read the cell probabilities (these will match the length of subset_barcodes)
  cell_probabilities <- h5_file[['droplet_latents/cell_probability']]$read()

  # 5. Close the HDF5 file connection
  h5_file$close()

# Create a tibble to hold the results
  cell_data <- tibble(
    Barcode = subset_barcodes,
    Cell_Probability = cell_probabilities
  )

  # Write the data to a tab-separated file using readr
  write_tsv(cell_data, output_txt_file)
