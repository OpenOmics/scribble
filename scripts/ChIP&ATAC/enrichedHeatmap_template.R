####################
#
# Name: enrichedHeatmap_template.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: March 18, 2024
# 
#
# Important Notes:
# This was written for a specific lab for a particular project so certain details are hard-coded,
#   but easily adaptable. They include:
#      1. Gene information comes from hg38
#      2. Inputs are gene symbols, not ensembl IDs
#      3. Plotting function is hard-coded for 4 samples
#      4. Focus is on TSS sites with equal distance on either side
#      5. All genes are reoriented so that the first nucleosome is to the right of the TSS
#
##############
# PACKAGES

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("EnrichedHeatmap", quietly = TRUE))
    BiocManager::install("EnrichedHeatmap")
if(!require(TxDb.Hsapiens.UCSC.hg38.knownGene))  
  BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
if(!require(org.Hs.eg.db))
  BiocManager::install('org.Hs.eg.db')

library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(org.Hs.eg.db)
library(EnrichedHeatmap)
library(circlize)

geneGR <- keepStandardChromosomes(genes(TxDb.Hsapiens.UCSC.hg38.knownGene),
             pruning.mode="coarse")

#############
# FUNCTIONS

getGeneCoords <- function(geneSymbols) {
   eid <- select(org.Hs.eg.db, geneSymbols, "ENTREZID", "SYMBOL")[["ENTREZID"]]
   geneInfo <- data.frame(geneSymbol=geneSymbols, EntrezID=eid)
   print(geneInfo)
   geneSubset <- geneGR[which(geneGR$gene_id %in% eid)]
   Idx <- sapply(geneSubset$gene_id, function(x) {grep(x, geneInfo$EntrezID)} )
   geneSubset$geneSymbol <- geneInfo$geneSymbol[Idx]
   return(geneSubset)
}

getGeneStarts <- function(geneCoords) {
   geneStarts <- geneCoords
   start(geneStarts) <- ifelse(strand(geneCoords) == '+',
                               start(geneCoords), end(geneCoords))
   end(geneStarts) <- start(geneStarts)
   return(geneStarts)
}

createMatrices <- function(bwFiles, inCoords, extendValue=10000) {
   wVal <- floor(extendValue *2 / 100)
   mats <- vector(mode="list", length=length(bwFiles))
   for (i in 1:length(mats)) {
     signal <- import(bwFiles[i])
     mats[[i]] <- normalizeToMatrix(signal, inCoords, w=wVal,
                      extend=c(extendValue,extendValue),
		      mean_mode="weighted", value_column="score")
   }
   rownames(mats[[length(mats)]]) <- inCoords$geneSymbol
   return(mats)
}

makePlot <- function(mats, col_fun, IDs, ymax) {
   EnrichedHeatmap(mats[[1]], column_title=IDs[1], col=col_fun,
                   name= "Avg RPGC",
                   top_annotation = HeatmapAnnotation(enriched =
		     anno_enriched(gp=gpar(col="navy"), ylim= c(0,ymax))) ) +
   EnrichedHeatmap(mats[[2]], column_title=IDs[2], col=col_fun,
                   show_heatmap_legend = FALSE,
                   top_annotation = HeatmapAnnotation(enriched =
		     anno_enriched(gp=gpar(col="navy"), ylim= c(0,ymax))) ) +
   EnrichedHeatmap(mats[[3]], column_title=IDs[3], col=col_fun,
                   show_heatmap_legend = FALSE,
                   top_annotation = HeatmapAnnotation(enriched =
		     anno_enriched(gp=gpar(col="navy"), ylim= c(0,ymax))) ) +
   EnrichedHeatmap(mats[[4]], column_title=IDs[4], col=col_fun,
                   show_heatmap_legend = FALSE, show_row_names=T,
                   top_annotation = HeatmapAnnotation(enriched =
		     anno_enriched(gp=gpar(col="navy"), ylim= c(0,ymax))) )
}

#############
# USAGE

bws <- list.files(pattern=".bw")

conditionIDs <- c("P1","P2","P3","P4")

geneSymbols <- c("TOX", "TBX21", "ZEB2", "IKZF3", "ZETB32", "NR4A2")

geneCoords <- getGeneCoords(geneSymbols)
geneStarts <- getGeneStarts(geneCoords)

mats <- createMatrices(bws, geneStarts, extendValue=10000)

ymax <- quantile(unlist(mats), 0.975)
col_fun = colorRamp2(quantile(unlist(mats), c(0.05, 0.6,0.85,0.925,0.975)),
      c("white", "aliceblue","lightskyblue1", "deepskyblue","navy"))

makePlot(mats, col_fun, conditionIDs, ymax)
