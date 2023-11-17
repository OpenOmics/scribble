# What is this

A utility script for downloading and verifying the integrity of single and paired end NGS reads from the SRA database. 
For a small number of SRR IDs this could be run through an interactive compute node serially, for a large collection of
SRR IDs this should be run through many small + quick compute jobs (swarm is the easiest option).

# How to run

How to run these scripts. These scripts are written specifically for the biowulf nih cluster and will need to be tailor to 
any other infrastructure that does not use lua environment modules to load software.

## get_srr_id.sh

__Before using__: set `data_dir` to your storage directory

This script is for single or serial collection of SRR ids. It will run `prefetch`, `vdb-validate`, and `fasterq-dump` in order then
compress the reads for storage.

Single sample:
```./get_srr_id.sh <SRRID>```

Example output:
```
--- downloading SRR24433054 ---
2023-11-17T22:15:58 prefetch.3.0.5: Current preference is set to retrieve SRA Normalized Format files with full base quality scores.
2023-11-17T22:15:58 prefetch.3.0.5: 1) Downloading '<SRRID>'...
2023-11-17T22:15:58 prefetch.3.0.5: SRA Normalized Format file is being retrieved, if this is different from your preference, it may be due to current file availability.
2023-11-17T22:15:58 prefetch.3.0.5:  Downloading via HTTPS...
2023-11-17T22:16:27 prefetch.3.0.5:  HTTPS download succeed
2023-11-17T22:16:29 prefetch.3.0.5:  '<SRRID>' is valid
...
```

## pull_sra.sh

```
./pull_sra.sh <SRRID1>,<SRRID2>,<SRRID3>,...,<SRRIDn>
```

Example output:
```
<swarm_job_id>
Swarm submitted
```

This script is for swarm/sbatch file creation for download and collection of large collection of SRR IDs, these file download concurrently.

# Integrity validation

`vdb-validate` output is collected and catalogued to ensure integrity of reads downloaded.
