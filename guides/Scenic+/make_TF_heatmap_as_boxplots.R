# This code was originally used in an Rmd script.
# The purpose is to take each Region AUC, Target Gene AUC, and TF expression and put them all together for each eRegulon as a single boxplot figure.

library(ggplot2)
library(patchwork)
# library(tidyr)

outfolder <- "scplus_CD4/outs/"

eRegulon <- read.delim(paste0(outfolder, "eRegulon_direct.tsv"))
eReg <- unique(eRegulon[,c("TF","Gene_signature_name","Region_signature_name","eRegulon_name", "regulation")])
rm(eRegulon)

region <- read.csv(paste0(outfolder, "eRegulon_region_AUC.csv"))
gene <- read.csv(paste0(outfolder, "eRegulon_gene_AUC.csv")
TF <- read.csv(paste0(outfolder, "eRegulon_TF_exp.csv")

names(region) <- c("Cell", eReg$eRegulon_name)
names(gene) <- c("Cell", eReg$eRegulon_name)
names(TF) <- c("Cell", eReg$eRegulon_name)

region2 <- tidyr::pivot_longer(region, !Cell, names_to="eRegulon", values_to="Value")
region2$Measurement <- "RegionAUC"
gene2	<- tidyr::pivot_longer(gene, !Cell, names_to="eRegulon", values_to="Value")
gene2$Measurement <- "TargetGeneAUC"
TF2	<- tidyr::pivot_longer(TF, !Cell, names_to="eRegulon", values_to="Value")
TF2$Measurement <- "TFexp"

allData<- rbind(region2, gene2, TF2)
allData$CellType <- gsub("_[0-9]+", "", allData$Cell)

for(i in 1:nrow(eReg)) {

tmp <-allData[which(allData$eRegulon == eReg$eRegulon_name[i]),]

p <- ggplot(tmp, aes(x=CellType, y=Value, color=CellType))
p <- p + #geom_violin() +
    geom_boxplot() +
    facet_wrap(. ~ Measurement, scales = "free") +
    theme_bw() +
    ggtitle(eReg$eRegulon_name[i])

print(p)
}
