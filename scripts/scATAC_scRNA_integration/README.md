scATAC & scRNA integration


integrate_scATAC_scRNA_liger.R
Integrate scATAC and scRNA using liger.

Samples were processed in Seurat/Signac using the basic parameters [Here's the signac vignette](https://stuartlab.org/signac/articles/pbmc_vignette) to get the "expression counts".

<br>
<br>
   
For running rliger on Biowulf as of 8/8/24:
`module load rliger`

This module is a singularity object with the most up to date versions of all required and suggested packages on the rliger cran website. 
It can be used interactively or in batch scripts, but cannot be used with Rstudio.

It has two run functions:  
`rliger-R` for running interactive commandline R  
`rliger-Rscript` for running R jobs in batch scripts  
