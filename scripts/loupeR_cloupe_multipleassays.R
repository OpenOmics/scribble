library(Seurat)
library(loupeR)
library(tidyverse)

# Script to create loupe browser files for ADT and GEX assays using loupeR
# in the new loupe browser version, you can view 2 objects side-by-side which is nice for visualizing ADT and GEX assays simultaneously (unable to do with cellranger aggr)

## Folder to save the objects to to
outdir <- 'loupe_export_folder'

## Reading in the Seurat object
seur <- readRDS('seur.rds')

# create loupe for GEX assay
create_loupe_from_seurat(seur, output_dir = outdir, 
                         output_name="loupeR_GEX")

# format metadata
metadata <- seur@meta.data %>%
  as.data.frame() 
metadata <- setNames(lapply(names(metadata), function(i) as.factor(metadata[[i]])), colnames(metadata))

# create loupe browser object with CITE-Seq assay using counts, umap projections, and metadata
create_loupe(
  seur@assays$ADT@counts,
  clusters = metadata,
  output_dir = outdir,
  projections = list(umap = seur@reductions$umap@cell.embeddings),
  output_name = 'loupeR_ADT'
)