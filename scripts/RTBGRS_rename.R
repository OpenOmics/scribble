Purpose: To semi-automate making useful filenames for RTBGRS
      data.

Written by: Tovah Markowitz
Date: 10/11/24

General idea: Takes a folder of fastq files from RTBGRS and
      their sample sheet as inputs. Assumes that the fastq
      files are named using the "Sample ID" column and
      names of interest are in "Sample Name". Makes an .sh
      script that when run will produces soft symlinks in
      the location where the script is run. Assumes all fastqs
      are gzipped.

Note: To prevent errors, inFolder should be a full path.

Example usage:
  inFolder <- "/data/RTB_GRS/IDSS_Projects/RTBGRS-105/fastq"
  inSampleSheet <- "GRS_0440_Samples_2024-09-27_14-51-01.xlsx"
  makeSymlinksScript(inFolder, inSampleSheet)



makeSymlinksScript <- function(inFolder, inSampleSheet) {
   library(readxl)
   inFiles <- list.files(inFolder, full.names = TRUE)
   # This is to determine if the file was R1 or R2
   R12 <- regmatches(inFiles, regexpr("R[12]+", inFiles) )
   # This is to extract the sample ID
   sampID <- gsub("_$", "", regmatches(inFiles, regexpr("LIB_[0-9_]*", inFiles) ) )
   info1 <- data.frame(inFiles, R12, sampID, gz)
 
   sampleInfo <- read_excel(inSampleSheet)
   sampleInfo2 <- data.frame(sampID=sampleInfo$`Sample ID`,
                             sampName=sampleInfo$`Sample Name`)


   info3 <- merge(info1, sampleInfo2)
   
   cat("#!/bin/bash", file="makeSymlinks.sh", sep="\n")
   for ( i in 1:nrow(info3) ) {
      tmp <- paste0("ln -s ", info3$inFiles[i], " ",
                    info3$sampName[i], ".", info3$R12[i],
                    ".fastq.gz\n")
      cat(tmp, file="makeSymlinks.sh", append=TRUE)
   }
}
