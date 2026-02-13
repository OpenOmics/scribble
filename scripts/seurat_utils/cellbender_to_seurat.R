#!/usr/bin/env Rscript

# CellBender H5 to Seurat v5 Converter
# This script reads a CellBender output .h5 file and creates a Seurat v5 object

# Load required libraries
library(Seurat)
library(hdf5r)
library(Matrix)

# Function to read CellBender H5 file
read_cellbender_h5 <- function(h5_file) {
  
  # Open the H5 file
  h5 <- H5File$new(h5_file, mode = "r")
  
  # List available groups to understand structure
  cat("H5 file structure:\n")
  print(names(h5))
  
  # CellBender typically stores data in 'matrix' group
  # Try different possible locations
  if (h5$exists("matrix")) {
    group <- h5[["matrix"]]
  } else {
    # List all groups to help debug
    cat("\nAvailable groups in file:\n")
    cat(paste(names(h5), collapse = "\n"), "\n")
    h5$close_all()
    stop("Could not find 'matrix' group in H5 file")
  }
  
  # Read the sparse matrix components
  data <- group[["data"]][1:group[["data"]]$dims]
  indices <- group[["indices"]][1:group[["indices"]]$dims]
  indptr <- group[["indptr"]][1:group[["indptr"]]$dims]
  shape <- group[["shape"]][1:group[["shape"]]$dims]
  
  # Read barcodes and features
  barcodes <- group[["barcodes"]][1:group[["barcodes"]]$dims]
  features <- group[["features"]][["name"]][1:group[["features"]][["name"]]$dims]
  
  # Close the H5 file
  h5$close_all()
  
  # Create sparse matrix (CellBender uses CSC format like 10X)
  mat <- sparseMatrix(
    i = indices + 1,  # R is 1-indexed
    p = indptr,
    x = as.numeric(data),
    dims = shape,
    index1 = TRUE
  )
  
  # Set column names (barcodes)
  colnames(mat) <- barcodes
  
  # Handle duplicate feature names
  if (any(duplicated(features))) {
    cat("Warning: Found", sum(duplicated(features)), "duplicate feature names\n")
    cat("Making feature names unique...\n")
    features <- make.unique(as.character(features))
  }
  
  # Set row names (features)
  rownames(mat) <- features
  
  return(mat)
}

# Function to convert CellBender H5 to Seurat object
cellbender_to_seurat <- function(h5_file, assay_name = "RNA", project_name = "CellBender") {
  
  # Check if file exists
  if (!file.exists(h5_file)) {
    stop("File not found: ", h5_file)
  }
  
  cat("Reading CellBender H5 file:", h5_file, "\n")
  
  # Read the H5 file using custom function
  counts <- read_cellbender_h5(h5_file)
  
  # Create Seurat object with v5 assay
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = project_name,
    assay = assay_name,
    min.cells = 0,  # Set filtering thresholds as needed
    min.features = 0
  )
  
  cat("Seurat object created successfully!\n")
  cat("Number of cells:", ncol(seurat_obj), "\n")
  cat("Number of features:", nrow(seurat_obj), "\n")
  
  return(seurat_obj)
}

# Main execution
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript cellbender_to_seurat.R <input.h5> [output.rds] [project_name]\n")
    cat("\nArguments:\n")
    cat("  input.h5      : Path to CellBender output H5 file\n")
    cat("  output.rds    : (Optional) Path to save Seurat object (default: input_seurat.rds)\n")
    cat("  project_name  : (Optional) Project name for Seurat object (default: CellBender)\n")
    quit(status = 1)
  }
  
  # Get arguments
  h5_file <- args[1]
  output_file <- if (length(args) >= 2) args[2] else gsub("\\.h5$", "_seurat.rds", h5_file)
  project_name <- if (length(args) >= 3) args[3] else "CellBender"
  
  # Convert to Seurat
  seurat_obj <- cellbender_to_seurat(h5_file, project_name = project_name)
  
  # Save the Seurat object
  cat("Saving Seurat object to:", output_file, "\n")
  saveRDS(seurat_obj, file = output_file)
  
  cat("\nConversion complete!\n")
  cat("Seurat object saved to:", output_file, "\n")
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}

# Example usage in R session:
# source("cellbender_to_seurat.R")
# seurat_obj <- cellbender_to_seurat("SCAF4785_corrected_filtered.h5")
# saveRDS(seurat_obj, "my_seurat_object.rds")
