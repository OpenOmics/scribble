# This script was written using R/4.5.0 on Biowulf.
# The purpose of this script is to make the equivalent of a DotPlot 
#    for a single gene when splitting/grouping on two variables, 
#    since the plot in Seurat requires each split.by value to be
#    associated with a unique color.
#
# The function called "BalloonPlot" requires 3 packages and
# has 4 input variables:
#      obj: the name of the Seurat object
#      gene: the name of the gene
#      group1: the name of the metadata column that will become the x-axis
#      group2: the name of the metadata column that will become the y-axis
#
# Note: No grouping values/metadata can contain underscores as the code is
#   currently written.

library(Seurat)
library(dplyr)
library(ggpubr)

BalloonPlot <- function(obj, gene, group1, group2) {
  avgExp <- AverageExpression(obj, features=gene,
     group.by=c(group1, group2))

  geneExp <- FetchData(obj, vars=c(gene, group1, group2))

  pctExp <- geneExp %>%
    group_by(.data[[group1]], .data[[group2]]) %>%
    summarise(
      pctExp = mean(.data[[gene]] != 0) * 100,
      .groups = 'drop'
    )

  names(pctExp)[1:2] <- c("group1", "group2")

  tmp <- colnames(avgExp$RNA)
  tmp <- matrix(unlist(strsplit(tmp,split="_")),ncol=2,byrow=T)
  plotDataA <- data.frame(group1=tmp[,1],
                       group2=tmp[,2],
                       avgExp=as.vector(avgExp$RNA))

  plotDataB <- merge(plotDataA, pctExp)

  print(ggballoonplot(plotDataB, x= "group1", y="group2",
                  size="pctExp", fill="avgExp",
                  size.range=c(0.1,20)) +
              scale_fill_viridis_c())
}
