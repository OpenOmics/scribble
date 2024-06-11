# General notes

As learned by Katie Hornick and Tovah Markowitz

NOTE: All scripts here are examples and need some tuning for each individual project, but notes are included in each to suggest where these changes need to take place.

## Updated notes as of 6/3/24:
1. The newest version of scenic+ now uses python/3.11 which is not available on Biowulf yet.
2. These documents are for the version that had been available back in 11/23 (now found in the github branch called "old")
3. The version created to use this documentation did not work with all samples analyzed. Reason unknown.
4. New versions of Ray do not allow more than 10 threads to be used simultaneously, does not like usage of lscratch, and can not have more than one job run in parallel.
5. To use the current version of scenic+, access https://github.com/OpenOmics/dockerfiles/tree/main/scenicplus/1.0a1.


## Tool documentation comes from:  
- https://pycistopic.readthedocs.io/en/latest/index.html  
- https://pycistarget.readthedocs.io/en/latest/index.html  
- https://scenicplus.readthedocs.io/en/latest/index.html 
  - https://github.com/aertslab/scenicplus/tree/old/notebooks

## Key information:  
The attached scripts take all the information from these various tool websites and streamlines them to one general pipeline.
This pipeline assumes the following:  
1. The dataset must be both scRNA and scATAC.  
    - If you only have scATAC, you can run all pycistopic and pycistarget steps, but none of the scenic+.  
    - If you only have scRNA, look into scenic instead of scenic+.  
2. This code is written to compare at least two conditions/cell types. If this is not true, you will need to adapt some steps.  
3. All conditions that are going to be compared should have the same labels and metadata header for scRNA and scATAC. You will only be able to look at one metadata column per analysis.  
4. This code assumes that the scATAC and scRNA are not from the same cells. If this is not true, you will need to adapt some steps in the scenic+ preparation steps.  


## Order of processing:  

1. Pull our [scenicplus docker image](https://hub.docker.com/r/skchronicles/scenicplusc) with singularity.
```bash
module load singularity
SINGULARITY_CACHEDIR=$PWD/.${USER} singularity pull -F docker://skchronicles/scenicplus:v0.1.0
```
2. seurat_scenic_prep.R: To convert seurat/signac objects to files that can be used by python  
3. pycistopic_model.py: part1 in processing scATAC data, typically takes 1-5 days to run and requires 200-500G of memory. Memory spike should happen within the first 2 hours. Each model will require 1 thread, but Ray only allows up to 10 threads. See details in python script. See submission_scripts.md for example run command.  
4. pycistopic_part2.py: part2 in processing of scATAC data. This step is very fast compared to the previous step, requires only 1 thread, and typically needs 50-200G of memory. See submission_scripts.md for example run command.  
5. a. pycistarget.py: part3 in processing if scATAC data only, requires 6 cpus, assume over an hour of run time. See submission_scripts.md for example run command. This is the last step if there is only scATAC data. Note: Not tested on dockerized version of the pipeline.  
   b. pycistarget_snakemake_prep.py: runs the first half of the pycistarget.py script to put pycistopic_part2.py outputs into the correct format for the snakemake pipeline  
6. Create custom cistarget database. This is a new step in the pipeline that we have not tested.  


6. scenic+_prep.py: combine information from pycistopic, pycistarget, and scRNA information into a single object  
7. Get the TF names and known target tables. See example below for human or try https://github.com/aertslab/scenicplus/tree/main/resources.  
`wget -O utoronto_human_tfs_v_1.01.txt  http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt`  
8. Set up bedToBigBed tool. Make sure its in a folder you can point to.  
`wget -O bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed`  
`chmod +x bedToBigBed`  
9. Run scenic+: This step only works on Biowulf when n_cpu=1 so it takes some time to run. It also uses a lot of memory; I ended up setting the lscratch to 900 for 20k scRNA and 20k scATAC cells and that still wasn't sufficient. If it needs more temp space, it will give a "Ray storage space" error. I ended up needing around 5T for 24h.
10. See https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html for good documentation on how to access the results, do some minor cleanup, and create some basic figures.  
