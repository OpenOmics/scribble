# Cell Ranger Aggregation Scripts

* `generateAggregateCSV_ATAC.py`
* `generateAggregateCSV_GEX.py`
* `generateAggregateCSV_MULTIOME.py`
* `generateAggregateCSV_VDJ.py`

These are basic scripts to help set up the csv file needed to run Cell Ranger aggregate. It should be run in the parent folder that contains all the `count` or `vdj` runs that will be aggregated. There is currently no script to assist with `multi` runs.

# Cell Ranger Summary Report Script

* `generateSummaryFiles.py`
* `generateSummaryFiles_ATAC.py`
* `generateSummaryFiles_MULTI.py`
* `generateSummaryFiles_MULTIOME.py`

These are scripts to help collect the statistics in the summary csv files for each sample into one Excel file, and also create a copy of each web summary file that has the sample name appended to the front. All of these files are gathered into a newly created finalreports folder.

`generateSummaryFiles.py` can work for GEX, VDJ, and feature barcode runs processed via `cellranger count`.

Currently the ATAC, MULTI, and MULTIOME version of the scripts try to filter out any aggregate runs by search for the intermediate folder of the cellranger `count` or `multi` run. That line would need to be commented out if those folders have been removed.
