# General notes

As learned by Katie Hornick and Tovah Markowitz

NOTE: All scripts here are examples and need some tuning for each individual project, but notes are included in each to suggest where these changes need to take place.

## Updated notes as of 6/3/24:
1. The newest version of scenic+ now uses python/3.11 which is not available on Biowulf yet.
4. New versions of Ray do not allow more than 10 threads to be used simultaneously and does not like usage of lscratch.
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
3. All conditions that are going to be compared should have the same labels and metadata header for scRNA and scATAC. You will only be able to look at one metadata column per analysis. In the example scripts, this is called "CellType".  
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
5. a. pycistarget.py: part3 in processing if scATAC data only, requires 6 cpus, assume over an hour of run time. See submission_scripts.md for example run command. This is the last step if there is only scATAC data. Note: Not tested on dockerized version of the pipeline. **OR**  
   b. pycistarget_snakemake_prep.py: runs the first half of the pycistarget.py script to put pycistopic_part2.py outputs into the correct format for the snakemake pipeline  
6. a. Create custom cistarget database. This is a new step in the pipeline that we have not tested. **OR**  
   b. Get required files from https://resources.aertslab.org/cistarget/. For hg38, they are:  
         - ctx_db file/rankings file: hg38_screen_v10_clust.regions_vs_motifs.rankings.feather  
         - dem_db file/scores file: hg38_screen_v10_clust.regions_vs_motifs.scores.feather  
         - motif_annotations: motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl  
7. Initialize the snakemake directory. See submission_scripts.md for example run command.  
8. Update config.yaml file and replace version in initialized directory in Snakemake/config subfolder. See config_notes.md for details on options.
9. Run snakemake. See submission_scripts.md for example run command.
10. To do downstream steps, follow website examples in an interactive session with the singularity object activated.
```bash
module load singularity
singularity exec -c -B $PWD scenicplus_v0.1.0.sif /bin/bash
```
