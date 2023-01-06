####################
#
# Name: convert_conns_into_cleaned-bedpe.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: January 6, 2023
# 
####################
#
# Purpose: To take a cicero conns object and convert it into a bedpe-like format
#          to use with the cicero flexdashboard I created.
#
# Note: will temporarily create a unique bedpe file around 200Gb in size
#
# Functions:
#    conns2Bedpe(conns, RDSfile)
#
# Variables:
#    conns		the cicero conns output data
#    RDSfile	the output name of the RDS to be created
#
####################


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("rtracklayer", quietly = TRUE))
    BiocManager::install("rtracklayer")

if (!require("InteractionSet", quietly = TRUE))
    BiocManager::install("InteractionSet")

library(InteractionSet)

conns2Bedpe <- function(conns, RDSfile) {
  tmp1 <- matrix(unlist(strsplit(conns$Peak1,"-")),ncol=3,byrow=T)
  tmp2 <- matrix(unlist(strsplit(as.character(conns$Peak2),"-")),ncol=3,byrow=T)
  conns2 <- data.frame(chrom1=tmp1[,1], start1=tmp1[,2], end1=tmp1[,3],
                     chrom2=tmp2[,1], start2=tmp2[,2], end2=tmp2[,3],
                     name=".", score=conns$coaccess)

  tmpfile <- tempfile(tmpdir=".",fileext=".bedpe")
  write.table(conns2,tmpfile,quote=F,row.names=F,sep="\t",col.names=F)

  conns <- rtracklayer::import(tmpfile)
  conns2 <- GInteractions(first(conns),second(conns))
  conns2$score <- mcols(conns)$score

  conns2 <- unique(swapAnchors(conns2))
  conns2 <- conns2[which(!is.na(conns2$score))]

  saveRDS(conns2, RDSfile)
  unlink(tmpfile)
}
