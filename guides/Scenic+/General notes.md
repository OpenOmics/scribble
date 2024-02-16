All scripts here are examples and need some tuning for each individual project, but notes are included in each to suggest where these changes need to take place.


Tool documentation comes from:  
https://pycistopic.readthedocs.io/en/latest/index.html  
https://pycistarget.readthedocs.io/en/latest/index.html  
https://scenicplus.readthedocs.io/en/latest/index.html


The attached scripts take all the information from these various tool websites and streamlines them to one general pipeline.
This pipeline assumes the following:  
1. The dataset must be both scRNA and scATAC.  
    - If you only have scATAC, you can run all pycistopic and pycistarget steps, but none of the scenic+.  
    - If you only have scRNA, look into scenic instead of scenic+.  
2. This code is written to compare at least two conditions/cell types. If this is not true, you will need to adapt some steps.  
3. All conditions that are going to be compared should have the same labels and metadata header for scRNA and scATAC. You will only be able to look at one metadata column per analysis.  
4. This code assumes that the scATAC and scRNA are not from the same cells. If this is not true, you will need to adapt some steps in the scenic+ preparation steps.  


Order of processing:  
1. prep_conda.sh: To set up the conda environment containing all necessary packages  
2. seurat_scenic_prep.R: To convert seurat/signac objects to files that can be used by python  
3. pycistopic_model.py: part1 in processing scATAC data, may need to be run multiple times, see note in script about changing cpu requirements, takes a number of hours to run
4. pycistopic_part2.py: part2 in processing of scATAC data
5. pycistarget.py: part3 in processing if scATAC data, requires 6 cpus, assume over an hour of run time
6. scenic+_prep.py: combine information from pycistopic, pycistarget, and scRNA information into a single object  
7. Get the TF names and known target tables. See example below for human or try https://github.com/aertslab/scenicplus/tree/main/resources.  
`wget -O utoronto_human_tfs_v_1.01.txt  http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt`  
8. Set up bedToBigBed tool  
`wget -O bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed`  
`chmod +x bedToBigBed`  
9. Run scenic+: This step only works on Biowulf when n_cpu=1 so it takes a lot of time to run. It also uses a lot of memory; I ended up setting the lscratch to 900 for 20k scRNA and 20k scATAC cells
