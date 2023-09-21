# Cell Ranger Aggregation 

Aggregation can only be performed if the same reference is used on all samples. If this is not the case, then the samples would need to be reprocessed to use the same reference. 

This is just a basic overview of how to set up aggregation runs. 

## Single Cell GEX (processed by `cellranger count`)

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

`cellranger aggr` has three arguments when running:
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
* sample_id: Identifier for the sample, this is what the generated reports and Loupe object will use as the name for the sample. It does not need to match the name used for when running `cellranger count`
* molecule_h5: Path to the molecule_info.h5 file produced by `cellranger count`

`generateAggregateCSV_GEX.py` is a very basic script that when run in a folder containing the `cellranger count` runs, it would find all the molecule_info.h5 files and generate a file named AggregatedDatasets.csv that can be used for aggregation

Additional information can be found at https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-3p-aggr

## Single Cell VDJ

### Example slurm job

```
#!/bin/sh
#SBATCH --job-name=countaggr_VDJ
#SBATCH --time=72:00:00
#SBATCH --mem=75g
#SBATCH --ntasks=18

module load cellranger

cellranger aggr --id=Aggregate --csv=AggregatedVDJ.csv --localcores=16 --localmem=32 2>aggr_run_v1.err 1>aggr_run_v1.log
```

### Information

`cellranger aggr` has two arguments when running:
* id - Required: Name of the resulting output
* csv - Required: CSV file containing the paths to the vdj_contig_info.pb files, donor, and origin information

The format of the CSV file is as follows
```
sample_id,vdj_contig_info,donor,origin
Spleen,/opt/runs/Spleen_vdj/outs/vdj_contig_info.pb,MM3,Spleen
AxLN,/opt/runs/AxLN_vdj/outs/vdj_contig_info.pb,MM3,AxLN
Liver,/opt/runs/Liver_vdj/outs/vdj_contig_info.pb,MM3,Liver
MeLN,/opt/runs/MeLN_vdj/outs/vdj_contig_info.pb,MM3,MeLN
PBMC,/opt/runs/PBMC_vdj/outs/vdj_contig_info.pb,MM3,PBMC
```
* sample_id: Identifier for the sample, this is what the generated reports and Loupe object will use as the name for the sample. It does not need to match the name used for when running `cellranger vdj`
* vdj_contig_info: Path to the vdj_contig_info.pb file produced by `cellranger vdj`
* donor: Identifier for the individual the sample was collected from
* origin: Specific source of the sample, which could be timepoint, tissue, or other metadata
* Optional other metadata columns can be added to describe the samples

`generateAggregateCSV_VDJ.py` is a very basic script that when run in a folder containing the `cellranger vdj` runs, it would find all the vdj_contig_info.pb files and generate a file named AggregatedDatasets.csv with the donor and origin column headers in place. This additional information will need to be filled in to use the file for aggregation

Additional information can be found at https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-5p-aggr

## Single Cell Feature Barcode or Hashtag (processed by `cellranger count` or `multi`)

### Information

The set up and required files are the same as those for the [single cell GEX run](#single-cell-gex-processed-by-cellranger-count). The additional thing to keep in mind is that aggregate will only work if the feature or CMO/hashtag reference is the same across all samples. This means that even if the sample does not contain those features or CMO/hashtags, it would need to be reprocessed while providing those as reference.

Additional information about how this could be run can be found at https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-3p-aggr#cmo-aggr

## Single Cell Multi with GEX and VDJ

### Example slurm job

```
#!/bin/sh
#SBATCH --job-name=countaggr_MULTI
#SBATCH --time=72:00:00
#SBATCH --mem=96g
#SBATCH --ntasks=18

module load cellranger

cellranger aggr --id=Aggregate --csv=AggregatedMulti.csv --localcores=16 --localmem=75 2>aggr_run_v1.err 1>aggr_run_v1.log
```

### Information

10x Genomics does not support the use of cellranger aggr to aggregate the outputs of T cell libraries enriched for gamma (TRG) and delta (TRD) chains.

`cellranger aggr` takes a CSV file that has a list of `cellranger multi` output directorys and can perform aggregation on different combinations of libraries (technologies), but it requires the same combination in all samples in order to be aggregated.

`cellranger aggr` has two arguments when running:
* id - Required: Name of the resulting output
* csv - Required: CSV file containing the paths to the per_sample_outs folder generated by `cellranger multi` 

The format of the CSV file is as follows
```
sample_id,sample_outs,donor,origin,VaccinationStatus
Sample1,/opt/runs/Sample1/outs/per_sample_outs/Sample1,D1,pbmc_t0,Pre-Vaccination
Sample2,/opt/runs/Sample2/outs/per_sample_outs/Sample2,D1,pbmc_t1,Post-Vaccination
```
* sample_id: Identifier for the sample, this is what the generated reports and Loupe object will use as the name for the sample. It does not need to match the name used for when running `cellranger multi`
* sample_outs: Path to the per_sample_outs folder produced by `cellranger multi`
* donor: Identifier for the individual the sample was collected from
* origin: Specific source of the sample, which could be timepoint, tissue, or other metadata
* Optional other metadata columns can be added to describe the samples

Additional information about this can be found at https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-5p-aggr#multi_aggr

## Single Cell ATAC

### Example slurm job

```
#!/bin/sh
#SBATCH --job-name=countaggr_ATAC
#SBATCH --time=72:00:00
#SBATCH --mem=75g
#SBATCH --ntasks=18

module load cellranger-atac

cellranger-atac aggr --id=Aggregate --csv=AggregatedATAC.csv --reference=/data/NCBR/references/cellranger_references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --localcores=16 --localmem=32 2>aggr_run_v1.err 1>aggr_run_v1.log
```

### Information

`cellranger-atac aggr` has three required arguments when running:
* id - Name of the resulting output
* csv - CSV file containing the paths to the fragments and single cell files.
* reference - Path to the folder containing the reference 

The format of the CSV file is as follows
```
library_id,fragments,cells
LV123,/opt/runs/LV123/outs/fragments.tsv.gz,/opt/runs/LV123/outs/singlecell.csv
LB456,/opt/runs/LB456/outs/fragments.tsv.gz,/opt/runs/LB456/outs/singlecell.csv
LP789,/opt/runs/LP789/outs/fragments.tsv.gz,/opt/runs/LP789/outs/singlecell.csv
```
* library_id: Identifier for the sample, this is what the generated reports and Loupe object will use as the name for the sample. It does not need to match the name used for when running `cellranger-atac count`
* fragments: Path to the fragments.tsv.gz file produced by `cellranger-atac count`
* cells: Path to the singlecell.csv file produced by `cellranger-atac count`

`generateAggregateCSV_ATAC.py` is a very basic script that when run in a folder containing the `cellranger-atac count` runs, it would find all the fragments.tsv.gz and singlecell.csv files and generate a file named AggregatedDatasets.csv that can be used for aggregation. The file would only contain the runs where it found both files.

Additional information can be found at https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/aggr

## Single Cell Multiome

### Example slurm job

```
#!/bin/sh
#SBATCH --job-name=countaggr_MULTIOME
#SBATCH --time=72:00:00
#SBATCH --mem=96g
#SBATCH --ntasks=18

module load cellranger-arc

cellranger-arc aggr --id=Aggregate --csv=AggregatedARC.csv --reference=/data/NCBR/references/cellranger_references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --localcores=16 --localmem=75 2>aggr_run_v1.err 1>aggr_run_v1.log
```

### Information

`cellranger-arc aggr` has three required arguments when running:
* id - Name of the resulting output
* csv - CSV file containing the paths to the fragments, barcodes, and molecule info files.
* reference - Path to the folder containing the reference 

The format of the CSV file is as follows
```
library_id,atac_fragments,per_barcode_metrics,gex_molecule_info
LV123,/opt/runs/LV123/outs/atac_fragments.tsv.gz,/opt/runs/LV123/outs/per_barcode_metrics.csv,/opt/runs/LV123/outs/gex_molecule_info.h5
LB456,/opt/runs/LB456/outs/atac_fragments.tsv.gz,/opt/runs/LB456/outs/per_barcode_metrics.csv,/opt/runs/LB456/outs/gex_molecule_info.h5
LP789,/opt/runs/LP789/outs/atac_fragments.tsv.gz,/opt/runs/LP789/outs/per_barcode_metrics.csv,/opt/runs/LP789/outs/gex_molecule_info.h5
```
* library_id: Identifier for the sample, this is what the generated reports and Loupe object will use as the name for the sample. It does not need to match the name used for when running `cellranger-arc count`
* atac_fragments: Path to the atac_fragments.tsv.gz file produced by `cellranger-arc count`
* per_barcode_metrics: Path to the per_barcode-metrics.csv file produced by `cellranger-arc count`
* gex_molecule_info: Path to the gex_molecule_info.h5 file produced by `cellranger-atac count`

`generateAggregateCSV_MULTIOME.py` is a very basic script that when run in a folder containing the `cellranger-arc count` runs, it would find all the atac_fragments.tsv.gz and identify all runs where the per_barcode_metrics.csv and gex_molecule_info.h5 files also exist. It will generate a file named AggregatedDatasets.csv containing the identified runs with all three files that can be used for aggregation. 

Additional information can be found at https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/aggr

