#!/usr/bin/env Rscript

# CellBender H5 to Seurat v5 Converter for Multiome Data
# This script reads a CellBender output .h5 file and creates a Seurat v5 object
# Supports multiome data (RNA + ATAC)

# Load required libraries
library(Seurat)
library(Signac)
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
  if (h5$exists("matrix")) {
    group <- h5[["matrix"]]
  } else {
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
  
  # Try to read feature types if available (for multiome)
  feature_types <- NULL
  if (group[["features"]]$exists("feature_type")) {
    feature_types <- group[["features"]][["feature_type"]][1:group[["features"]][["feature_type"]]$dims]
    cat("Found feature types:\n")
    print(table(feature_types))
  }
  
  # Close the H5 file
  h5$close_all()
  
  # Create sparse matrix (CellBender uses CSC format like 10X)
  mat <- sparseMatrix(
    i = indices + 1,
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
  
  return(list(matrix = mat, feature_types = feature_types))
}

# Function to convert CellBender H5 to Seurat object (multiome)
cellbender_to_seurat_multiome <- function(h5_file, 
                                          atac_h5 = NULL,
                                          fragments = NULL,
                                          genome = NULL,
                                          project_name = "Multiome") {
  
  if (!file.exists(h5_file)) {
    stop("File not found: ", h5_file)
  }
  
  cat("Reading CellBender H5 file (RNA):", h5_file, "\n")
  
  # Read the CellBender corrected RNA data
  cb_data <- read_cellbender_h5(h5_file)
  rna_counts <- cb_data$matrix
  feature_types <- cb_data$feature_types
  
  # If feature types exist, separate RNA and ATAC
  if (!is.null(feature_types)) {
    cat("\nSeparating RNA and ATAC modalities...\n")
    
    rna_idx <- feature_types == "Gene Expression"
    atac_idx <- feature_types == "Peaks"
    
    # Subset the original matrix for both modalities
    full_counts <- rna_counts
    rna_counts <- full_counts[rna_idx, ]
    atac_counts <- full_counts[atac_idx, ]
    
    cat("RNA features:", nrow(rna_counts), "\n")
    cat("ATAC features:", nrow(atac_counts), "\n")
    
    # Create Seurat object with RNA
    seurat_obj <- CreateSeuratObject(
      counts = rna_counts,
      project = project_name,
      assay = "RNA",
      min.cells = 0,
      min.features = 0
    )
    
    # Add ATAC assay if peaks exist
    if (nrow(atac_counts) > 0) {
      cat("Adding ATAC assay...\n")
      
      if (!is.null(fragments) || !is.null(genome)) {
        seurat_obj[["ATAC"]] <- CreateChromatinAssay(
          counts = atac_counts,
          sep = c(":", "-"),
          fragments = fragments,
          genome = genome
        )
      } else {
        cat("Creating ATAC assay without fragments (can be added later)\n")
        seurat_obj[["ATAC"]] <- CreateAssayObject(counts = atac_counts)
      }
    }
    
  } else {
    # No feature types found - assume RNA only from CellBender
    cat("\nNo feature types found - treating as RNA only\n")
    
    seurat_obj <- CreateSeuratObject(
      counts = rna_counts,
      project = project_name,
      assay = "RNA",
      min.cells = 0,
      min.features = 0
    )
    
    # Add ATAC from separate file if provided
    if (!is.null(atac_h5)) {
      cat("Reading ATAC data from:", atac_h5, "\n")
      atac_counts <- Read10X_h5(atac_h5)
      
      if (is.list(atac_counts)) {
        atac_counts <- atac_counts$Peaks
      }
      
      if (!is.null(fragments) || !is.null(genome)) {
        seurat_obj[["ATAC"]] <- CreateChromatinAssay(
          counts = atac_counts,
          sep = c(":", "-"),
          fragments = fragments,
          genome = genome
        )
      } else {
        cat("Creating ATAC assay without fragments (can be added later)\n")
        seurat_obj[["ATAC"]] <- CreateAssayObject(counts = atac_counts)
      }
    }
  }
  
  cat("\nSeurat object created successfully!\n")
  cat("Number of cells:", ncol(seurat_obj), "\n")
  cat("Assays:", names(seurat_obj@assays), "\n")
  
  return(seurat_obj)
}

# Main execution
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript cellbender_to_seurat.R <rna_cellbender.h5> [atac.h5] [fragments.tsv.gz] [genome] [output.rds] [project_name]\n")
    cat("\nArguments:\n")
    cat("  rna_cellbender.h5 : Path to CellBender corrected RNA H5 file\n")
    cat("  atac.h5           : (Optional) Path to ATAC H5 file if separate\n")
    cat("  fragments.tsv.gz  : (Optional) Path to fragments file for ATAC\n")
    cat("  genome            : (Optional) Genome annotation (e.g., 'hg38', 'mm10')\n")
    cat("  output.rds        : (Optional) Path to save Seurat object\n")
    cat("  project_name      : (Optional) Project name for Seurat object\n")
    quit(status = 1)
  }
  
  h5_file <- args[1]
  atac_h5 <- if (length(args) >= 2 && args[2] != "NA") args[2] else NULL
  fragments <- if (length(args) >= 3 && args[3] != "NA") args[3] else NULL
  genome <- if (length(args) >= 4 && args[4] != "NA") args[4] else NULL
  output_file <- if (length(args) >= 5) args[5] else gsub("\\.h5$", "_seurat.rds", h5_file)
  project_name <- if (length(args) >= 6) args[6] else "Multiome"
  
  seurat_obj <- cellbender_to_seurat_multiome(
    h5_file, 
    atac_h5 = atac_h5,
    fragments = fragments,
    genome = genome,
    project_name = project_name
  )
  
  cat("Saving Seurat object to:", output_file, "\n")
  saveRDS(seurat_obj, file = output_file)
  
  cat("\nConversion complete!\n")
  cat("Seurat object saved to:", output_file, "\n")
}

if (!interactive()) {
  main()
}

# Example usage:
# source("cellbender_to_seurat.R")
# seurat_obj <- cellbender_to_seurat_multiome("SCAF4785_MEN_2077_corrected_filtered.h5")
# saveRDS(seurat_obj, "multiome_seurat.rds")
