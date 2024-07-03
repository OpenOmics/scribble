####################
#
# Name: Heatmap_aroundTSS.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: December 6, 2023
#
####################
#
# Example code to:
# 1. Extract data from specific windows from multiple bigwigs.
#     In this case, we are extracting the information at TSS sites
#     and ordering all sites to have the first codon to the right of
#     of the TSS.
# 2. Create a heatmap of the results.


library(GenomicRanges)
library(EnrichedHeatmap)
library(circlize)
library(rtracklayer)

bws <- list.files(path="bigwig",pattern=".Q5DD.RPGC.inputnorm.bw",full.names=T)
bwNames <- c("A","B","C","D","E","F")

# inBed needs to be a GenomicRanges object of just genes going into the heatmap
bed_starts <- inBed
start(bed_starts) <- ifelse(strand(inBed) == '+', start(inBed), end(inBed))
end(bed_starts) <- start(bed_starts)

mats <- vector(mode="list",length=length(bws)) 

for (i in 1:length(bws)) {
  signal <- import(bws[i])

  mats[[i]] <- normalizeToMatrix(signal, bed_starts,
          extend=c(10000, 10000), w=250,
          # adjust extend and w as needed where extend is the distance upstream and 
          # downstream of the granges window in order, and w is the window size to be
          # averaged. Use a higher window size to decrease the size of the final image
          # as files can become large but greater than 300 bins per inch on a graphic
          # tends to be useless
          mean_mode="weighted", value_column="score")
}

rownames(mats[[5]]) <- bed_starts$name

# this needs to be tuned to the data
col_fun = colorRamp2(quantile(unlist(mats), c(0.05, 0.5,0.8,0.9,0.95)), c("white", "aliceblue","lightskyblue1", "deepskyblue","navy"))

# this will automatically make the heatmap and the profile plots
# ylim is associated with the profile plot
# name is the legend label
# sigGenes2 is a vector of the labels of the genes in the order seen in inBed
#  For example: "significant" and "not significant"
# col=2:3 has two values as there are two values in sigGenes2
EnrichedHeatmap(mats[[1]], row_split=sigGenes2,column_title=bwNames[1],col = col_fun, name="K4Me3",
    top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 2:3), ylim=c(-5,100)))) +
EnrichedHeatmap(mats[[2]], row_split=sigGenes2,column_title=bwNames[2],col = col_fun, show_heatmap_legend = FALSE,
    top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 2:3), ylim=c(-5,100))) ) +
EnrichedHeatmap(mats[[3]], row_split=sigGenes2,column_title=bwNames[3],col = col_fun, show_heatmap_legend = FALSE,
    top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 2:3), ylim=c(-5,100))) ) +
EnrichedHeatmap(mats[[4]], row_split=sigGenes2,column_title=bwNames[4],col = col_fun, show_heatmap_legend = FALSE,
    top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 2:3), ylim=c(-5,100))) ) +
EnrichedHeatmap(mats[[5]], row_split=sigGenes2,column_title=bwNames[5],col = col_fun, show_heatmap_legend = FALSE,
    top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 2:3), ylim=c(-5,100))) ) +
EnrichedHeatmap(mats[[6]], row_split=sigGenes2,column_title=bwNames[6],col = col_fun, show_heatmap_legend = FALSE,
    top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 2:3), ylim=c(-5,100))), show_row_names=T)
