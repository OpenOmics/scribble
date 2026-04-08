This code was shared with me by Cihan and involves code on both python and R.

Python:

import scanpy as sc
import pandas as pd
from scipy import io

adata = sc.read_h5ad("your_data.h5ad")
# Export sparse matrix
io.mmwrite("counts.mtx", adata.X.T) 
# Export metadata and gene names
adata.obs.to_csv("metadata.csv")
pd.DataFrame(adata.var_names).to_csv("genes.csv", index=False)



R:

library(Seurat)
library(Matrix) 
counts <- readMM("counts.mtx")
metadata <- read.csv("metadata.csv", row.names = 1)
genes <- read.csv("genes.csv") 
rownames(counts) <- genes[,1]
colnames(counts) <- rownames(metadata) 
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
