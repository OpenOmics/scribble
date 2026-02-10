#!/usr/bin/env Rscript

# Script to run SoupX on GEX portion of 10X Multiome data - Version 2
# Required packages: SoupX, Seurat, hdf5r, optparse

suppressPackageStartupMessages({
  library(optparse)
  library(SoupX)
  library(Seurat)
  library(hdf5r)
})

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

option_list <- list(
  make_option(c("-f", "--filtered"),
              type = "character",
              default = NULL,
              help = "Path to filtered multiome .h5 file [required]",
              metavar = "FILE"),
  
  make_option(c("-r", "--raw"),
              type = "character",
              default = NULL,
              help = "Path to raw feature-barcode matrix .h5 file [required]",
              metavar = "FILE"),
  
  make_option(c("-o", "--output"),
              type = "character",
              default = "soupx_output",
              help = "Output directory [default: %default]",
              metavar = "DIR"),
  
  make_option(c("-c", "--contamination"),
              type = "double",
              default = NULL,
              help = "Manual contamination fraction (0-1). If not specified, will auto-estimate",
              metavar = "FLOAT"),
  
  make_option(c("-d", "--dims"),
              type = "integer",
              default = 10,
              help = "Number of PCA dimensions for clustering [default: %default]",
              metavar = "INT"),
  
  make_option(c("-s", "--resolution"),
              type = "double",
              default = 0.8,
              help = "Clustering resolution [default: %default]",
              metavar = "FLOAT"),
  
  make_option(c("-m", "--markers"),
              type = "character",
              default = "CD3D,CD79A",
              help = "Comma-separated marker genes for diagnostic plots [default: %default]",
              metavar = "GENES"),
  
  make_option(c("-p", "--prefix"),
              type = "character",
              default = "soupX",
              help = "Prefix for output files [default: %default]",
              metavar = "STRING"),
  
  make_option(c("-v", "--verbose"),
              action = "store_true",
              default = FALSE,
              help = "Print verbose output")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nRun SoupX ambient RNA correction on 10X Multiome GEX data",
  epilogue = paste(
    "\nExample usage:",
    "  Rscript soupx_multiome.R -f filtered.h5 -r raw.h5 -o output_dir",
    "  Rscript soupx_multiome.R -f filtered.h5 -r raw.h5 -c 0.1 -m CD3E,MS4A1",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$filtered) || is.null(opt$raw)) {
  print_help(opt_parser)
  stop("Both --filtered and --raw arguments are required", call. = FALSE)
}

if (!file.exists(opt$filtered)) {
  stop(paste("Filtered h5 file not found:", opt$filtered), call. = FALSE)
}

if (!file.exists(opt$raw)) {
  stop(paste("Raw h5 file not found:", opt$raw), call. = FALSE)
}

if (!is.null(opt$contamination)) {
  if (opt$contamination < 0 || opt$contamination > 1) {
    stop("Contamination fraction must be between 0 and 1", call. = FALSE)
  }
}

# Parse marker genes
marker_genes <- strsplit(opt$markers, ",")[[1]]
marker_genes <- trimws(marker_genes)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

log_msg <- function(msg) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"))
}

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

log_msg("=== Starting SoupX analysis - Version 2 ===")
log_msg(paste("Filtered h5:", opt$filtered))
log_msg(paste("Raw h5:", opt$raw))
log_msg(paste("Output directory:", opt$output))

# Create output directory
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

# Read the filtered (cells) multiome .h5 file and extract GEX
log_msg("Reading filtered multiome .h5 file...")
gex_data <- Read10X_h5(opt$filtered, use.names = TRUE)

# If the file contains multiple modalities, extract just GEX
if (is.list(gex_data)) {
  log_msg("Multiple modalities detected. Extracting Gene Expression data...")
  toc <- gex_data[["Gene Expression"]]
} else {
  toc <- gex_data
}
log_msg(paste("Filtered matrix:", nrow(toc), "genes x", ncol(toc), "cells"))

# Read raw counts (all droplets including empty)
log_msg("Reading raw feature-barcode matrix...")
tod <- Read10X_h5(opt$raw, use.names = TRUE)

if (is.list(tod)) {
  tod <- tod[["Gene Expression"]]
}
log_msg(paste("Raw matrix:", nrow(tod), "genes x", ncol(tod), "droplets"))

# Create SoupChannel object
log_msg("Creating SoupChannel object...")
sc <- SoupChannel(tod, toc)

# Estimate contamination
# Perform quick clustering using Seurat
log_msg("Performing clustering for contamination estimation...")
log_msg(paste("  Using", opt$dims, "PCA dimensions"))
log_msg(paste("  Clustering resolution:", opt$resolution))

seurat_obj <- CreateSeuratObject(counts = toc)
seurat_obj <- NormalizeData(seurat_obj, verbose = opt$verbose)
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = opt$verbose)
seurat_obj <- ScaleData(seurat_obj, verbose = opt$verbose)
seurat_obj <- RunPCA(seurat_obj, verbose = opt$verbose)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:opt$dims, verbose = opt$verbose)
seurat_obj <- FindClusters(seurat_obj, resolution = opt$resolution, verbose = opt$verbose)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:opt$dims, verbose = opt$verbose)

log_msg(paste("Identified", length(unique(seurat_obj$seurat_clusters)), "clusters"))

# Add clustering info to SoupChannel
sc <- setClusters(sc, setNames(seurat_obj$seurat_clusters, colnames(seurat_obj)))

# Add dimension reduction for visualization
sc <- setDR(sc, seurat_obj@reductions$umap@cell.embeddings)

# Estimate contamination fraction
if (is.null(opt$contamination)) {
  log_msg("Auto-estimating contamination fraction...")
  sc <- autoEstCont(sc)
  global_rho <- mean(sc$metaData$rho)
  log_msg(paste("Estimated contamination:", round(global_rho * 100, 2), "%"))
} else {
  log_msg(paste("Using manual contamination fraction:", opt$contamination))
  sc <- setContaminationFraction(sc, opt$contamination)
  global_rho <- opt$contamination
}

# Calculate corrected counts
log_msg("Calculating corrected counts...")
out <- adjustCounts(sc)

# Save results
log_msg("Saving results...")

# Save corrected count matrix
saveRDS(out, file.path(opt$output, paste0(opt$prefix, "_corrected_counts.rds")))

# Create Seurat object with corrected counts
seurat_corrected <- CreateSeuratObject(counts = out, project = "SoupX_corrected")
saveRDS(seurat_corrected, file.path(opt$output, paste0(opt$prefix, "_seurat_corrected.rds")))

# Save contamination information - UPDATED FORMAT
log_msg("Writing contamination information...")

# Create structured contamination info
contamination_df <- data.frame(
  metric = c("global_rho", "n_cells", "n_genes", "n_clusters"),
  value = c(
    global_rho,
    ncol(sc$toc),
    nrow(sc$toc),
    length(unique(sc$metaData$clusters))
  ),
  stringsAsFactors = FALSE
)

write.csv(contamination_df, 
          file.path(opt$output, paste0(opt$prefix, "_contamination_info.csv")),
          row.names = FALSE)

log_msg(paste("Saved contamination info CSV with", nrow(contamination_df), "rows"))

# Save detailed summary - CLEANED UP VERSION
log_msg("Writing contamination summary...")
sink(file.path(opt$output, paste0(opt$prefix, "_contamination_summary.txt")))
cat("SoupX Contamination Estimation Summary\n")
cat("========================================\n\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste("Global contamination fraction (rho):", round(mean(sc$metaData$rho) * 100, 2), "%\n"))
cat(paste("Number of cells:", ncol(sc$toc), "\n"))
cat(paste("Number of genes:", nrow(sc$toc), "\n"))
cat(paste("Number of clusters:", length(unique(sc$metaData$clusters)), "\n\n"))

cat("Per-cluster contamination:\n")
cluster_rho <- tapply(sc$metaData$rho, sc$metaData$clusters, mean)
# Check if contamination varies by cluster
if (length(unique(cluster_rho)) == 1) {
  cat("  Note: Global contamination applied uniformly across all clusters\n")
  cat(sprintf("  All clusters: %.3f%%\n\n", cluster_rho[1] * 100))
} else {
  for (i in seq_along(cluster_rho)) {
    cat(sprintf("  Cluster %s: %.3f%%\n", names(cluster_rho)[i], cluster_rho[i] * 100))
  }
  cat("\n")
}

# Extract top soup genes more robustly
if (!is.null(sc$soupProfile)) {
  cat("Top 10 genes in soup profile:\n")
  # Handle different soupProfile structures
  if (is.matrix(sc$soupProfile) || is.data.frame(sc$soupProfile)) {
    soup_vals <- sc$soupProfile[, 1]
    names(soup_vals) <- rownames(sc$soupProfile)
  } else if (is.vector(sc$soupProfile)) {
    soup_vals <- sc$soupProfile
  } else {
    soup_vals <- NULL
  }
  
  if (!is.null(soup_vals) && length(soup_vals) > 0) {
    top_soup <- head(sort(soup_vals, decreasing = TRUE), 10)
    for (i in seq_along(top_soup)) {
      cat(sprintf("  %d. %s: %.4f\n", i, names(top_soup)[i], top_soup[i]))
    }
  } else {
    cat("  (Soup profile structure not recognized)\n")
  }
} else {
  cat("Top 10 genes in soup profile:\n")
  cat("  (No soup profile available)\n")
}
sink()

log_msg("Summary file written")

# Generate gene-specific change maps
log_msg("Generating gene-specific change maps...")

# Find top expressed genes as backup markers
top_genes <- names(sort(Matrix::rowSums(sc$toc), decreasing = TRUE)[1:20])

# Check which requested markers are available
available_markers <- marker_genes[marker_genes %in% rownames(sc$toc)]
if (length(available_markers) == 0) {
  log_msg("Warning: None of the specified marker genes found in data")
  log_msg(paste("Using top 3 expressed genes instead:", paste(head(top_genes, 3), collapse = ", ")))
  plot_genes <- head(top_genes, 3)
} else {
  log_msg(paste("Found", length(available_markers), "of", length(marker_genes), "specified markers"))
  plot_genes <- available_markers
}

# Generate change map plots - using manual method for reliability
plots_generated <- 0
pdf_path <- file.path(opt$output, paste0(opt$prefix, "_change_maps.pdf"))

pdf(pdf_path, width = 14, height = 6)

for (gene in plot_genes) {
  tryCatch({
    log_msg(paste("  Plotting change map for", gene, "..."))
    
    if (!gene %in% rownames(sc$toc)) {
      log_msg(paste("    Gene", gene, "not in matrix, skipping"))
      next
    }
    
    # Get expression before and after
    expr_before <- as.numeric(sc$toc[gene, ])
    expr_after <- as.numeric(out[gene, ])
    
    # Get UMAP coordinates from Seurat object
    umap_coords <- seurat_obj@reductions$umap@cell.embeddings
    
    # Ensure we have matching dimensions
    if (nrow(umap_coords) != length(expr_before)) {
      log_msg(paste("    Dimension mismatch: UMAP has", nrow(umap_coords), "cells but expression has", length(expr_before)))
      next
    }
    
    # Handle zeros and create better color scaling
    expr_before_nonzero <- expr_before[expr_before > 0]
    expr_after_nonzero <- expr_after[expr_after > 0]
    
    # Create color palette
    colors <- colorRampPalette(c("grey90", "lightblue", "blue", "red", "darkred"))(100)
    
    # Create side-by-side plots
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
    
    # Before plot
    if (length(expr_before_nonzero) > 0) {
      max_expr <- max(expr_before)
      breaks <- seq(0, max_expr, length.out = 101)
      color_idx <- cut(expr_before, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      color_idx[is.na(color_idx)] <- 1
    } else {
      color_idx <- rep(1, length(expr_before))
    }
    
    plot(x = umap_coords[, 1], 
         y = umap_coords[, 2],
         col = colors[color_idx],
         pch = 16, cex = 0.5,
         xlab = "UMAP_1", ylab = "UMAP_2",
         main = paste0(gene, " - Before SoupX\n(mean: ", round(mean(expr_before), 2), ")"))
    
    # After plot
    if (length(expr_after_nonzero) > 0) {
      max_expr <- max(expr_after)
      breaks <- seq(0, max_expr, length.out = 101)
      color_idx <- cut(expr_after, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      color_idx[is.na(color_idx)] <- 1
    } else {
      color_idx <- rep(1, length(expr_after))
    }
    
    plot(x = umap_coords[, 1], 
         y = umap_coords[, 2],
         col = colors[color_idx],
         pch = 16, cex = 0.5,
         xlab = "UMAP_1", ylab = "UMAP_2",
         main = paste0(gene, " - After SoupX\n(mean: ", round(mean(expr_after), 2), ")"))
    
    plots_generated <- plots_generated + 1
    log_msg(paste("    SUCCESS"))
    
  }, error = function(e) {
    log_msg(paste("    FAILED:", e$message))
  })
}

dev.off()

# Check if PDF was actually created with content
if (plots_generated == 0) {
  log_msg("Warning: No change maps were generated")
  if (file.exists(pdf_path)) {
    file.remove(pdf_path)
    log_msg("Removed empty PDF file")
  }
} else {
  pdf_size <- file.info(pdf_path)$size
  log_msg(paste("Generated", plots_generated, "change map(s) in", basename(pdf_path)))
  log_msg(paste("PDF size:", format(pdf_size, big.mark = ","), "bytes"))
  
  # More reasonable threshold for PDFs with actual plots
  if (pdf_size < 10000) {
    log_msg("Warning: PDF file seems small, please verify it opens correctly")
  }
}

# Save session info
writeLines(capture.output(sessionInfo()), 
           file.path(opt$output, paste0(opt$prefix, "_session_info.txt")))

log_msg("=== SoupX Analysis Complete ===")
log_msg(paste("Contamination fraction:", round(global_rho * 100, 2), "%"))
log_msg(paste("Results saved to:", opt$output))
log_msg("\nOutput files:")
log_msg(paste("  -", paste0(opt$prefix, "_corrected_counts.rds: Corrected count matrix")))
log_msg(paste("  -", paste0(opt$prefix, "_seurat_corrected.rds: Seurat object with corrected counts")))
log_msg(paste("  -", paste0(opt$prefix, "_contamination_info.csv: Contamination statistics")))
log_msg(paste("  -", paste0(opt$prefix, "_contamination_summary.txt: Detailed contamination summary")))
if (plots_generated > 0) {
  log_msg(paste("  -", paste0(opt$prefix, "_change_maps.pdf: Gene-specific change maps (", plots_generated, " genes)")))
}
log_msg(paste("  -", paste0(opt$prefix, "_session_info.txt: R session information")))
