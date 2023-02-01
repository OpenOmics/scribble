####################
#
# Name: regioneR_bootstrap.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: November 30, 2022
#
####################
#
# Purpose: To make regioneR easier to use with functions to help load the
#          necessary data and prep it for analysis, and an example for the most
#          common usage cases
#
# See bottom of this script for example usage.
#


####################
# Load required packages

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
     install.packages("BiocManager")
    }
  BiocManager::install("GenomicRanges")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
     install.packages("BiocManager")
    }
    BiocManager::install("rtracklayer")
}

if (!requireNamespace("regioneR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("regioneR")
}


####################
# Functions for prepping files for bootstrapping

expandRegions <- function(grObject, finalSize) {
  # This function will resize all rows in a GenomicRanges object around their
  # midpoints to reach a desired final size.
  library(GenomicRanges)
  midpoint <- floor(width(grObject) / 2)
  start(grObject) <- start(grObject) + midpoint
  end(grObject) <- start(grObject)
  grObject <- shift(grObject, shift=-finalSize/2)
  grObject <- resize(grObject, width=finalSize)
  return(grObject)
}

prepFile <- function(inFile, addFlanks=0) {
  # This function will import a bed or gff file and prepare it for downstream
  # analyses. Steps include importing the data as a GenomicRanges object,
  # sorting it, extending the regions as requested, and collapsing overlapping
  # rows.
  library(GenomicRanges)
  fileData <- rtracklayer::import(inFile)
  fileData <- sort(fileData)
  if (addFlanks != 0) {
    fileData <- regioneR::extendRegions(fileData, extend.start=addFlanks, extend.end=addFlanks)
  }
  return(fileData)
}

###################
# Specialized functions for evaluation

overlapAmt <- function(A,B,...) {
  sum(width(commonRegions(A,B)))
}

###################
# Example usage

data1 <- prepFile("file1.bed", addFlanks=(1e5/2))
data2 <- expandRegions(prepFile("file2.bed"), finalSize=1000)
blacklist <- prepFile("hg19.blacklist.bed")
chrs <- read.table("hg19.len", header=FALSE)
# chrs should be a two column file of chromosome name and chromosome length



library(regioneR)
pt <- permTest(A=data1, B=data2, ntimes=1000, genome=chrs, mask=blacklist,
               # randomization options
               randomize.function=randomizeRegions, allow.overlaps=FALSE,
               # evaluation options
               evaluate.function=numOverlaps, count.once=TRUE,
               # statistical test option
               alternative="auto")
plot(pt)

# Options for evaluation.function:
    # numOverlaps: number of overlapping features
        # count.once: each feature in A is counted only once
    # meanDistance: mean distance between features in A nearest feature in B
    # overlapAmt: number of bases overlapping between the two sets of features
