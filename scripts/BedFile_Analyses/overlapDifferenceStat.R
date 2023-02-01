####################
#
# Name: overlapDifferenceStat.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: February 1, 2023
#
####################
#
# Purpose: Have two peak files and want to see if one of them
#          overlaps a feature type of interest (like promoters)
#          more often than the other? And determine if the
#          difference is significant? That's what this script is
#          for. You can also add flanks to the features of interest
#          to see if one of the files is more likely to be within
#          a certain distance of the features of interest.
#
# Inputs:  All input files must work with rtracklayer, meaning
#          that they should be in bed, gtf, or gff format.
#
# Outputs: It will print the contingency table and p-value to
#          screen.
#          
# Function:
#    overlapDifferenceStat(file1, file2, overlapFile, addFlanks=0)
#
# Variables:
#    file1:       Name of the first file. 
#    file2:       Name of the second file.
#    overlapFile: Name of the file with the features to be overlapped
#                 against. In other words, we will measure the overlap
#                 the number of times a feature in this file is
#                 overlapped by features in each of the other two files.
#    addFlanks:   To make the features of the overlapFile wider. This value
#                 in base pairs will be added to both sides of the features
#                 prior to doing the overlap analysis.
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
# Functions 

library(regioneR)

prepFile <- function(inFile, addFlanks=0) {
  # This function will import a bed or gff file and prepare it for downstream
  # analyses. Steps include importing the data as a GenomicRanges object,
  # sorting it, and extending the regions as requested.
  fileData <- rtracklayer::import(inFile)
  fileData <- sort(fileData)
  if (addFlanks != 0) {
    fileData <- extendRegions(fileData, extend.start=addFlanks, extend.end=addFlanks)
  }
  return(fileData)
}

makeContingencyTable <- function(data1, data2, data3) {
  ov1 <- overlapRegions(data3, data1, only.boolean = TRUE)
  ov2 <- overlapRegions(data3, data2, only.boolean = TRUE)
  overlap1=factor(ov1,levels=c(TRUE,FALSE))
  overlap2=factor(ov2,levels=c(TRUE,FALSE))
  mat1 <- table(data.frame(overlap1,overlap2))
  return(mat1)
}

overlapDifferenceStat <- function(file1, file2, overlapFile, addFlanks=0) {
  data1 <- prepFile(file1)
  data2 <- prepFile(file2)
  data3 <- prepFile(overlapFile, addFlanks=addFlanks)
  cTable <- makeContingencyTable(data1, data2, data3)
  print(cTable)
  tmp <- fisher.test(cTable)
  print(paste0("The p-value using the fisher exact test is: ", tmp$p.value ))
}
