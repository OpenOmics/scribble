Note: All input_data and output_data should be absolute paths.  

Only mentioning parameters that have obvious reasons to change. See website for all other parameters.  

## input_data ##  
**cisTopic_obj_fname:**  path to the pickle object produced by pycistopic_part2.py on line 38 (pycistopic.pkl)  
**GEX_anndata_fname:** path to the scRNA h5ad file  
**region_set_folder:** path to the bedfiles folder produced by pycistarget_snakemake_prep.py  
**ctx_db_fname:** path to the rankings file (.feather)  
**dem_db_fname:** path to the scores file (.feather)  
**path_to_motif_annotations:** path to the motif annotations file (.tbl)  

## output_data ##  
Make sure to change the absolute path to go into the outs directory in the snakemake folder or it will end up in the Snakemake subfolder  

## params_general ##  
**temp_dir:**  
**n_cpu:** keep this under the number of cpus requested in the snakemake sbatch or steps will error out. I used 8 here and 20 in the sbatch.  

## parames_data_preparation ##  
**is_multiome:** set to True or False  
**key_to_group_by:** If not multiome, change this to the variable name used to compare along. I used "CellType" in the examples.  
**nr_cell_per_metacells:** Typically kept at 10 for 5-20k cells in analysis. Lower to 5 if very few cells; raise if more to lower memory consumption if problematic  
**species:** "hsapiens" or "mmusculus" for example  
**biomart_host:** link to the version of ensembl that best matches the reference genome used in your project  

## params_motif_enrichment ##  
**species:** "homo_sapiens" or "mus_musculus" for example  
**annotation_version:** change this depending on the name of the motif annotations file  
