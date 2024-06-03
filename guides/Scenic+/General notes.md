# General notes

As learned by Katie Hornick and Tovah Markowitz

All scripts here are examples and need some tuning for each individual project, but notes are included in each to suggest where these changes need to take place.

**Updated notes as of 6/3/24:**
1. The newest version of scenic+ now uses python/3.11 which is not available on Biowulf yet.
2. These documents are for the version that had been available back in 11/23 (now found in the github branch called "old")
3. The version created to use this documentation did not work with all samples analyzed. Reason unknown.
4. New versions of Ray do not allow more than 10 threads to be used simultaneously, does not like usage of lscratch, and can not have more than one job run in parallel.
5. To use the current version of scenic+, access https://github.com/OpenOmics/dockerfiles/tree/main/scenicplus/1.0a1.


**Tool documentation comes from:** 
- https://pycistopic.readthedocs.io/en/latest/index.html  
- https://pycistarget.readthedocs.io/en/latest/index.html  
- https://scenicplus.readthedocs.io/en/latest/index.html 
  - https://github.com/aertslab/scenicplus/tree/old/notebooks


The attached scripts take all the information from these various tool websites and streamlines them to one general pipeline.
This pipeline assumes the following:  
1. The dataset must be both scRNA and scATAC.  
    - If you only have scATAC, you can run all pycistopic and pycistarget steps, but none of the scenic+.  
    - If you only have scRNA, look into scenic instead of scenic+.  
2. This code is written to compare at least two conditions/cell types. If this is not true, you will need to adapt some steps.  
3. All conditions that are going to be compared should have the same labels and metadata header for scRNA and scATAC. You will only be able to look at one metadata column per analysis.  
4. This code assumes that the scATAC and scRNA are not from the same cells. If this is not true, you will need to adapt some steps in the scenic+ preparation steps.  


## Getting started

**Order of processing:**  

1. Pull our [scenicplus docker image](https://hub.docker.com/r/skchronicles/scenicplusc) with singularity.
```bash
module load singularity
SINGULARITY_CACHEDIR=$PWD/.${USER} singularity pull -F docker://skchronicles/scenicplus:v0.1.0
```
2. seurat_scenic_prep.R: To convert seurat/signac objects to files that can be used by python  
3. pycistopic_model.py: part1 in processing scATAC data, may need to be run multiple times, see note in script about changing cpu requirements, takes a number of hours/days to run.

Here is an example of how to setup and run this on Biowulf:
```bash
# Create a tmp directory for this run, 
# so ray has its own dedicated workspace
mkdir -p tmp;
export tmp="$(mktemp -d -p "$PWD/tmp")";
echo "Created tmp directory: $tmp"

# Create SLURM jobs script to run 
# the pycistopic_model_CD8.py script
cat << EOF > run_scenicplus_models.sh
#!/usr/bin/env bash
#SBATCH --job-name=scenicplus
#SBATCH --mail-type=END,FAIL
#SBATCH --time=2-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10

set -e
module load singularity;

# Change to the directory where the 
# pycistopic_model_CD8.py script is
# located. The input files to the 
# script should also be in the same 
# directory as the script.
cd /data/path/to/project/scenic/;

echo "Starting to run scenicplus script"
# Please do not increase the number of
# CPUs in the pycistopic_model_CD8.py 
# script to more than 8. Increasing the 
# number of CPUs will cause ray workers 
# to unexpectedly error out. There is 
# an unresolved bug in ray that causes
# this issue. 
singularity exec -c -B $PWD,${tmp}:/tmp scenicplus_v0.1.0.sif /bin/bash -c "cd $PWD; python $PWD/pycistopic_model_CD8.py"
echo "Exit-code of scenicplus: $?"
EOF

# Submit the job to the cluster
chmod +x run_scenicplus_models.sh
sbatch run_scenicplus_models.sh
```

4. pycistopic_part2.py: part2 in processing of scATAC data

Here is an example of how to setup and run this on Biowulf:
```bash
# Create a tmp directory for this run, 
# so ray has its own dedicated workspace
mkdir -p tmp;
export tmp="$(mktemp -d -p "$PWD/tmp")";
echo "Created tmp directory: $tmp"

# Create SLURM jobs script to run 
# the pycistopic_part2.py script
cat << EOF > run_scenicplus_part2.sh
#!/usr/bin/env bash
#SBATCH --job-name=scenicplus
#SBATCH --mail-type=END,FAIL
#SBATCH --time=2-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=8

set -e
module load singularity;

# Change to the directory where the 
# pycistopic_part2.py script is
# located. The input files to the 
# script should also be in the same 
# directory as the script.
cd /data/path/to/project/scenic/;

echo "Starting to run scenicplus script"
singularity exec -c -B $PWD,${tmp}:/tmp scenicplus_v0.1.0.sif /bin/bash -c "cd $PWD; python $PWD/pycistopic_part2.py"
echo "Exit-code of scenicplus: $?"
EOF

# Submit the job to the cluster
chmod +x run_scenicplus_part2.sh
sbatch run_scenicplus_part2.sh
```

5. pycistarget.py: part3 in processing if scATAC data, requires 6 cpus, assume over an hour of run time

Here is an example of how to setup and run this on Biowulf:
```bash
# Create a tmp directory for this run, 
# so ray has its own dedicated workspace
mkdir -p tmp;
export tmp="$(mktemp -d -p "$PWD/tmp")";
echo "Created tmp directory: $tmp"

# Create SLURM jobs script to run 
# the pycistarget.py script
cat << EOF > run_scenicplus_pycistarget.sh
#!/usr/bin/env bash
#SBATCH --job-name=scenicplus
#SBATCH --mail-type=END,FAIL
#SBATCH --time=2-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=8

set -e
module load singularity;

# Change to the directory where the 
# pycistarget.py script is located. 
# The input files to the script 
# should also be in the same 
# directory as the script.
cd /data/path/to/project/scenic/;

echo "Starting to run scenicplus script"
# Please do not increase the number of
# CPUs in the pycistarget.py 
# script to more than 8. Increasing the 
# number of CPUs will cause ray workers 
# to unexpectedly error out. There is 
# an unresolved bug in ray that causes
# this issue. 
singularity exec -c -B $PWD,${tmp}:/tmp scenicplus_v0.1.0.sif /bin/bash -c "cd $PWD; python $PWD/pycistarget.py"
echo "Exit-code of scenicplus: $?"
EOF

# Submit the job to the cluster
chmod +x run_scenicplus_pycistarget.sh
sbatch run_scenicplus_pycistarget.sh
```
6. scenic+_prep.py: combine information from pycistopic, pycistarget, and scRNA information into a single object  
7. Get the TF names and known target tables. See example below for human or try https://github.com/aertslab/scenicplus/tree/main/resources.  
`wget -O utoronto_human_tfs_v_1.01.txt  http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt`  
8. Set up bedToBigBed tool. Make sure its in a folder you can point to.  
`wget -O bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed`  
`chmod +x bedToBigBed`  
9. Run scenic+: This step only works on Biowulf when n_cpu=1 so it takes some time to run. It also uses a lot of memory; I ended up setting the lscratch to 900 for 20k scRNA and 20k scATAC cells and that still wasn't sufficient. If it needs more temp space, it will give a "Ray storage space" error. I ended up needing around 5T for 24h.
10. See https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html for good documentation on how to access the results, do some minor cleanup, and create some basic figures.  
