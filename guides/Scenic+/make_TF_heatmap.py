# This needs to be run before you can use make_TF_heatmap.R
# Not all folder locations are accurate
# scenicplus_v0.1.0.sif needs to be active for some parts of this script to work

import mudata
scplus_mdata = mudata.read("scplusmdata.h5mu")

scplus_mudata = scplus_mdata

size_features, color_features, feature_names = scplus_mudata.uns['direct_e_regulon_metadata'][
            ["Region_signature_name", "Gene_signature_name", "eRegulon_name"]] \
            .drop_duplicates().values.T
# grabbing unique names

size_matrix = scplus_mudata["direct_region_based_AUC"].to_df()
color_matrix = scplus_mudata["direct_gene_based_AUC"].to_df()
group_by = scplus_mudata.obs["scRNA_counts:CCR"].tolist()

size_matrix = size_matrix[size_features]
color_matrix = color_matrix[color_features]

color_matrix.to_csv("eRegulon_gene_AUC.csv")
size_matrix.to_csv("eRegulon_region_AUC.csv")

# Calculate mean by group_by variable
color_matrix_avg = color_matrix.groupby(group_by).mean()
size_matrix_avg = size_matrix.groupby(group_by).mean()

color_matrix_avg.to_csv("eRegulon_gene_mean_AUC.csv")
size_matrix_avg.to_csv("eRegulon_region_mean_AUC.csv")

# This is to get the TF expression
from scenicplus.scenicplus_class import mudata_to_scenicplus
scplus_obj = mudata_to_scenicplus(
     mdata = scplus_mdata,
     path_to_cistarget_h5 = "ctx_results.hdf5",
     path_to_dem_h5 = "dem_results.hdf5")

color_matrix = scplus_obj.to_df('EXP')

size_matrix_features = (
        size_matrix.columns.to_list(),                      #full eRegulon name
        [f.split('_(')[0] for f in size_matrix.columns],    #eRegulon name without number of targets
        [f.split('_')[0] for f in size_matrix.columns])     #tf name

color_matrix = color_matrix[size_matrix_features[2]]

color_matrix.to_csv("eRegulon_TF_exp.csv")
