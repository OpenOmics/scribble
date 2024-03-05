library(Seurat)
library(dplyr)


## Script is assuming that the column Sample in the seurat object has Cell Ranger sample name
aggr <- read.csv("AggregateDatasets/outs/aggregation.csv")

## Folder to save the metadata to
outdir <- 'aggr_seurat_export'

## Reading in the Seurat object
seur <- readRDS('seur.rds')

## Listing metadata columns to keep. In this situation it is keeping the Azimuth predicted labels and the Seurat generated clusters
keep_cols <- c('Sample', grep('score', colnames(seur@meta.data)[c(grep('cr.', colnames(seur@meta.data)), grep('celltype', colnames(seur@meta.data)))], invert=TRUE, value=TRUE), grep('res.', colnames(seur@meta.data), value=TRUE))

## Extract data to export
metadata <- seur@meta.data[,keep_cols]
umap <- seur@reductions$umap@cell.embeddings

## Rename rows based on aggregation order
rows <- rownames(umap)
## Looking for which row in the aggr.csv file matches with the Sample name for that specific cell barcode. Renames the cell barcode index based on the row number. This is how the barcodes are renamed for Seurat integration. If the barcode is renamed in a different way, such as using new.cell.ids, this would need to be edited.
rows <- sapply(rows, function(x) (x %>% gsub(pattern='1', replacement=which(aggr$sample_id == seur$Sample[x])) %>% strsplit(split='_'))[[1]][1])
  

## Rename rows in data 
rownames(metadata) <- rows
rownames(umap) <- rows

## Put data into correct format for loupe browser
umap = cbind(rownames(umap), umap)
colnames(umap) = c('Barcode', 'UMAP-1', 'UMAP-2')
metadata = cbind('barcode' = rownames(metadata), metadata)

write.table(umap, file = file.path(outdir, paste0(patient, '_merge_umap.csv')), quote = F, sep = ',', row.names = F, col.names = T)
write.table(metadata, file = file.path(outdir, paste0(patient, '_metadata.csv')), quote = F, sep = ',', row.names = F, col.names = T)
