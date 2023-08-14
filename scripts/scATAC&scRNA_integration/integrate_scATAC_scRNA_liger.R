library(Signac)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(rliger)
## combined here is the scATAC-seq object
# scATAC-Seq data was processed om Signac using basic parameters
gene.activities <- GeneActivity(combined)
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)
DefaultAssay(combined) <- 'RNA'

mice.rna=readRDS("mice.rna.rds")
DefaultAssay(mice.rna) <- 'RNA'
data <- list(atac = as.matrix(combined@assays$RNA@counts), rna = as.matrix(mice.rna@assays$RNA@counts))


int <- createLiger(data)
int <- rliger::normalize(int)
int <- selectGenes(int)
int <- scaleNotCenter(int)
int <- optimizeALS(int, k = 20)
int <- quantile_norm(int)
int <- louvainCluster(int, resolution = 0.3,verbose = TRUE)
#int <- louvainCluster(int, resolution = 1,verbose = TRUE)
int <- runUMAP(int, distance = 'cosine', n_neighbors = 10, min_dist = 0.3)
save(int,file="liger_integration_umap.RData")

plots <- plotByDatasetAndCluster(int, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
plots[[1]] + plots[[2]]
png("liger_umaps.png", height=1300, width=3500, res=300)
plots[[1]] + plots[[2]]
dev.off()

## convert liger back to seurat object
int2 <- ligerToSeurat(int,by.dataset=T)
table(int2@active.ident)

int2$liger_clusters <- int2@active.ident
png("liger_split.png", height=1300, width=3500, res=300)
DimPlot(int2,group.by="liger_clusters",label=T,repel=T,split.by="orig.ident") + NoLegend()
dev.off()
