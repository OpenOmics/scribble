These are the basic steps for scATAC analyses as of March 21, 2023:  

Look at https://github.com/OpenOmics/scribble/blob/main/scripts/scATAC/scATAC_functions.R for example code for everything downstream of clustering.

1. When only working with one sample, follow https://stuartlab.org/signac/articles/pbmc_vignette.html or
   https://stuartlab.org/signac/articles/mouse_brain_vignette.html for initial processing until cluster identification
   for human or mouse, respectively.
2. If you have an scRNA-seq dataset you want to use for cluster annotations, also use the above html vignettes to "Create a gene activity matrix"
   and "Integrate with scRNA-seq". This is a basic label-transfer method. Please note that gene activity matrices should not be
   treated as if they are actual/relative gene expression patterns because they are not that accurate.
3. If you have multiple samples, follow the main instructions on https://stuartlab.org/signac/articles/merging.html for sample merging. Do not use
   section called "Merging without a common feature set" unless you have no other choice. If integration is needed, follow the same steps until the
   actual merge at which point you can integrate and do downstream steps.
4. After identifying clusters (and potentially annotating them), identify differential peaks.
5. For annotations, set up a txdb of only protein-coding genes and use ChIPseeker 
   https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html to make annotations and create visualizations.
6. Create a volcano plot of the data, using EnhancedVolcano or plotly. Optional.
7. Identify motifs enriched in the differential peaks.
