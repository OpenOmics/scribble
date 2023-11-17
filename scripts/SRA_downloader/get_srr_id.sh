#!/bin/bash
set -euo pipefail

function err() { cat <<< "$@" 1>&2; }
function fatal() { cat <<< "$@" 1>&2; exit 1; }
function retry() {
    # Tries to run a cmd 5 times before failing
    # If a command is successful, it will break out of attempt loop
    # Failed attempts are padding with the following exponential
    # back-off strategy {4, 16, 64, 256, 1024} in seconds
    # @INPUTS "$@"" = cmd to run
    # @CALLS fatal() if command cannot be run in 5 attempts
    local n=1
    local max=5
    local attempt=true # flag for while loop
    while $attempt; do
        # Attempt command and break out of attempt loop if successful
        "$@" && attempt=false || {
            # Try again up to 5 times
            if [[ $n -le $max ]]; then
                err "Command failed: $@"
                delay=$(( 4**$n ))
                err "Attempt: ${n}/${max}. Trying again in ${delay} seconds!\n"
                sleep $delay;
                ((n++))
            else
                fatal "Fatal: command failed after max retry attempts!"
            fi
        }
    done
}

SRR_ID=$1
data_dir=$(pwd)
lscratch=/lscratch/$SLURM_JOB_ID
module purge
module load sratoolkit
echo "--- downloading $SRR_ID ---"
#retry fasterq-dump --split-files --split-3 -F fastq -O /lscratch/$SLURM_JOB_ID --temp /lscratch/$SLURM_JOB_ID/tmp -p --threads 32 $SRR_ID
cd $lscratch
retry prefetch $SRR_ID
vdb-validate $SRR_ID >> $data_dir/md5s/${SRR_ID}.validate
fasterq-dump --split-files --split-3 -F fastq -p --threads 16 $SRR_ID
cd $data_dir
if [ -f /lscratch/$SLURM_JOB_ID/${SRR_ID}_1.fastq ]; then
    echo "/lscratch/$SLURM_JOB_ID/${SRR_ID}_1.fastq" > /data/NIAMS_IDSS/dev/NIAMS-40/${SRR_ID}.md5
    md5sum /lscratch/$SLURM_JOB_ID/${SRR_ID}_1.fastq >> /data/NIAMS_IDSS/dev/NIAMS-40/${SRR_ID}.md5
    echo >> /data/NIAMS_IDSS/dev/NIAMS-40/${SRR_ID}.md5
    echo "/lscratch/$SLURM_JOB_ID/${SRR_ID}_2.fastq" >> /data/NIAMS_IDSS/dev/NIAMS-40/${SRR_ID}.md5
    md5sum /lscratch/$SLURM_JOB_ID/${SRR_ID}_2.fastq >> /data/NIAMS_IDSS/dev/NIAMS-40/${SRR_ID}.md5
    echo >> /data/NIAMS_IDSS/dev/NIAMS-40/${SRR_ID}.md5
    sed -i -r "s/(^[\@\+]SRR\S+)/\1\/1/" /lscratch/$SLURM_JOB_ID/${SRR_ID}_1.fastq
    sed -i -r "s/(^[\@\+]SRR\S+)/\1\/1/" /lscratch/$SLURM_JOB_ID/${SRR_ID}_2.fastq
    pigz -c -9 -p 32 /lscratch/$SLURM_JOB_ID/${SRR_ID}_1.fastq > /data/NIAMS_IDSS/dev/NIAMS-40/SRA/${SRR_ID}_1.fastq.gz
    pigz -c -9 -p 32 /lscratch/$SLURM_JOB_ID/${SRR_ID}_2.fastq > /data/NIAMS_IDSS/dev/NIAMS-40/SRA/${SRR_ID}_2.fastq.gz
else
    echo "/lscratch/$SLURM_JOB_ID/${SRR_ID}.fastq" > /data/NIAMS_IDSS/dev/NIAMS-40/${SRR_ID}.md5
    md5sum /lscratch/$SLURM_JOB_ID/${SRR_ID}.fastq >> /data/NIAMS_IDSS/dev/NIAMS-40/md5s/${SRR_ID}.md5
    sed -i -r "s/(^[\@\+]SRR\S+)/\1\/1/" /lscratch/$SLURM_JOB_ID/${SRR_ID}.fastq
    pigz -c -9 -p 10 /lscratch/$SLURM_JOB_ID/${SRR_ID}.fastq > /data/NIAMS_IDSS/dev/NIAMS-40/SRA/${SRR_ID}.fastq.g
fi
rm /lscratch/$SLURM_JOB_ID/*.fastq
