#!/bin/bash
set -euo pipefail

for srr in $@; do
    echo "./get_srr_id.sh ${srr}"
done >> download_from_sra.sh

chmod +x download_from_sra.sh

# Submit the SRA fastq download
# script as a job to cluster
swarm --gres=lscratch:500 -b 25 -g 8 -t 16 download_from_sra.sh
echo "Swarm submitted"
