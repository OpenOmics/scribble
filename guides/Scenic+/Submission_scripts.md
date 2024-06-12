This is a list of the submission scripts needed to run the pipeline as found in: docker://skchronicles/scenicplus:v0.1.0

Scripts are in order of usage and indicate suggested sbatch submission guidelines, especially for Biowulf. They are not all steps in the workflow. 
See General Notes.md for the list of all steps.

## Run pycistopic_model.py

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
#SBATCH --time=5-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10

set -e
module load singularity;

# Change to the directory where the 
# pycistopic_model.py script is
# located. The input files to the 
# script should also be in the same 
# directory as the script.
cd /data/path/to/project/scenic/;

echo "Starting to run scenicplus script"
# Please do not increase the number of
# CPUs in the pycistopic_model.py 
# script to more than 8. Increasing the 
# number of CPUs will cause ray workers 
# to unexpectedly error out. There is 
# an unresolved bug in ray that causes
# this issue. 
singularity exec -c -B $PWD,${tmp}:/tmp scenicplus_v0.1.0.sif /bin/bash -c "cd $PWD; python $PWD/pycistopic_model.py"
echo "Exit-code of scenicplus: $?"
EOF

# Submit the job to the cluster
chmod +x run_scenicplus_models.sh
sbatch run_scenicplus_models.sh
```

## Run pycistopic_part2.py

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
#SBATCH --time=6:00:00
#SBATCH --mem=200G

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

## Run pycistarget.py
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

## Run pycistarget_snakemake_prep.py
```bash

# Create script to run 
# the pycistarget_snakemake_prep.py script
cat << EOF > run_scenicplus_pycistarget_snakemake_prep.sh
#!/bin/bash
set -e
module load singularity;

# Change to the directory where the 
# pycistarget_snakemake_prep.py script is located. 
# The input files to the script 
# should also be in the same 
# directory as the script.
cd /data/path/to/project/scenic/;

# make sure this matches hardcorded directory name in pycistarget_snakemake_prep.py
mkdir bedfiles

echo "Starting to run scenicplus script"
singularity exec -c -B $PWD,${tmp}:/tmp scenicplus_v0.1.0.sif /bin/bash -c "cd $PWD; python $PWD/pycistarget_snakemake_prep.py"
echo "Exit-code of scenicplus: $?"
EOF

# run as part of an interactive job, needs less than 1G ram and 1 thread
chmod +x run_scenicplus_pycistarget.sh
sh run_scenicplus_pycistarget_snakemake_prep.sh
```

## Initialize snakemake directory
```bash

# Create script to initialize
cat << EOF > initialize_snakemake.sh
#!/bin/bash

ml singularity

# Change to the directory where
# the input files are located.
cd /data/path/to/project/scenic/;

function scenicplus() { singularity exec -c -B $PWD,$PWD/tmp:/tmp scenicplus_v0.1.0.sif /usr/local/bin/scenicplus "$@" ; }
export -f scenicplus

DIR="scplus_snakemake"

mkdir $DIR
cd $DIR
mkdir tmp
mkdir outs

cd ..
scenicplus init_snakemake --out_dir $PWD/$DIR
EOF

# run as part of an interactive job, needs less than 1G ram and 1 thread
chmod +x initialize_snakemake.sh
sh initialize_snakemake.sh
```

