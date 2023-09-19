# Cell Ranger Aggregation 

## Single Cell GEX (processed by Cell Ranger count)

### Example slurm job

```
#!/bin/sh
#SBATCH --job-name=countaggr
#SBATCH --time=72:00:00
#SBATCH --mem=75g
#SBATCH --ntasks=18

module load cellranger

cellranger aggr --id=Aggregate --csv=AggregatedDatasets.csv --normalize=none --localcores=16 --localmem=32 2>aggr_run_v1.err 1>aggr_run_v1.log
```

### Information

Cell Ranger aggr has three arguments when running:
* id - Required: Name of the resulting output
* csv - Required: CSV file containing the paths to the molecule_info.h5 files.
* normalize - Optional: String specifying how to normalize depths across the different libraries, with two valid options - mapped (default) and none. 


The format of the CSV file is as follows
```
sample_id,molecule_h5
LV123,/opt/runs/LV123/outs/molecule_info.h5
LB456,/opt/runs/LB456/outs/molecule_info.h5
LP789,/opt/runs/LP789/outs/molecule_info.h5
```
* sample_id: Identifier for the sample, this is what the generated reports and Loupe object will use as the name for the sample
* molecule_h5: Path to the molecule_info.h5 file produced by Cell Ranger count

generateAggregateCSV.py is a very basic script that when run in a folder containing the cellranger count runs, it would find all the molecule_info.h5 files and generate a file named AggregatedDatasets.csv that can be used for aggregation

Additional information can be found at https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-3p-aggr
