# BedFile_Analyses

#### regioneR_bootstrap.R
> An R script designed to make regioneR (https://bioconductor.org/packages/release/bioc/html/regioneR.html) easier to use.
regioneR offers a statistical framework based on customizable permutation tests to assess the association between genomic 
region sets and other genomic features. This script includes code to help load the two region sets into R and prepare them
for analysis. It also includes an example for comparing the number of overlaps between the two region sets and notes on
a few other common comparisons that can be used instead.

#### CreatingIGVbatchFromBedFiles.py
> To create an IGV batch script from a bed file or bed-like file that will allow an semi-automated creation of snapshots across a genome. Just run
this script, open IGV, load the proper tracks and format them as you chose and then run the batch script to get a folder of pngs named by snapshop region.

