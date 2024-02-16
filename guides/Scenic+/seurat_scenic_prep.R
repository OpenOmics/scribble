# Name: seurat_scenic_prep.R
# Author: Tovah Markowitz
# Purpose: To take seurat/signac objects and convert them into files that will be useful for scenic+ analyses
# Note: these functions were written for Seurat version 4 and may need to be updated for future versions.

#####################################
# Usage

## Step 1. Load scATAC object
## Step 2. Subset to only include cells of interest
## Step 3. ct2 <- convertSeuratATAC(ATAC)
##         write.csv(ct2, "peak_counts.csv", quote=F)
##         ATACm <- data.frame(ATAC@meta.data)
##         write.table(ATACm,"ATAC_metadata.txt",quote=F,sep="\t")
## Step 4. Load scRNA object
## Step 5. Subset to only include cells of interest
## Step 6. convertSeuratRNA(RNA,"scRNA")
##         RNAm <- data.frame(RNA@meta.data)
##         write.table(RNAm,"RNA_metadata.txt",quote=F,sep="\t")
## Step 7: Optional: only required when data has been integrated
##               Need to have counts, data, and scale.data or information can end up 
##               with the incorrect label in the new format.
##               see: https://mojaveazure.github.io/seurat-disk/reference/Convert.html
##         convertSeuratRNA2(RNA,"scRNA2")

#####################################
# Packages

library(Seurat)
library(Signac)
library(SeuratDisk)

######################################
# Functions

convertSeuratATAC <- function(scATAC) {
  ct <- scATAC@assays$ATAC@counts
  dense_cistopic_count_matrix <- as.data.frame(as(ct, "matrix"))
  rownames(dense_cistopic_count_matrix) <- sub("-", ":", rownames(dense_cistopic_count_matrix), fixed = TRUE)
  return(dense_cistopic_count_matrix)
}

convertSeuratRNA <- function(rna, outroot) {
  i <- sapply(rna@meta.data, is.factor)
  rna@meta.data[i] <- lapply(rna@meta.data[i], as.character)

  SaveH5Seurat(rna, filename = paste0(outroot,".h5Seurat"), overwrite = T)
  Convert(paste0(outroot,".h5Seurat"), dest = "h5ad", overwrite = T,assay = DefaultAssay(object =rna))
}

convertSeuratRNA2 <- function(rna, outroot) {
  rna2 <- DietSeurat(rna)
  i <- sapply(rna2@meta.data, is.factor)
  rna2@meta.data[i] <- lapply(rna2@meta.data[i], as.character)

  SaveH5Seurat(rna2, filename = paste0(outroot,".h5Seurat"), overwrite = T)
  Convert(paste0(outroot,".h5Seurat"), dest = "h5ad", overwrite = T,assay = DefaultAssay(object =rna2))
}
