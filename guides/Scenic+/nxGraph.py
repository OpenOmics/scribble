### For this script to run and produce a figure, you must:
#          1) Create a subfolder to work in.
#          2) Place nxTable.py within the subfolder along with your version of this script.
#          3) Run the following in the outermost folder that contains all scenic+ data/codes: 
#                  module load singularity; singularity exec -c -B $PWD scenicplus_v0.1.0.sif /bin/bash
#          4) Move to the subfolder and run your version of this script.

# Note: This code is designed to start with the scenic+ version 2 object and then revert to the scenic+ version 1 
# object for plotting purposes. Code in nxTable.py is from the scenic+ version 1 code with minor adaptations to deal with
# key metadata changes between the two versions. See for more details: https://github.com/aertslab/scenicplus/issues/393

import mudata
scplus_mdata = mudata.read("../scplus_IL22IFNG/outs/scplusmdata.h5mu")

from scenicplus.scenicplus_class import mudata_to_scenicplus

scplus_obj = mudata_to_scenicplus(
    mdata = scplus_mdata,
    path_to_cistarget_h5 = "../scplus_IL22IFNG/outs/ctx_results.hdf5",
    path_to_dem_h5 = "../scplus_IL22IFNG/outs/dem_results.hdf5"
)

# needs filtering to remove the "extended eregulons"
direct_md = scplus_obj.uns['eRegulon_metadata']
direct_filt = direct_md[direct_md['Gene_signature_name'].str.contains("direct")]
scplus_obj.uns['direct_e_regulon_metadata_filtered'] = direct_filt

from pycisTopic.diff_features import find_highly_variable_features

hvr = find_highly_variable_features(scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata']['Region']))], n_top_features=1000, plot = False)
hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata']['Gene']))].T, n_top_features=1000, plot = False)

from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx, export_to_cytoscape
from nxTable import *
import matplotlib.pyplot as plt
import os

nx_tables = create_nx_tables(
    scplus_obj = scplus_obj,
    eRegulon_metadata_key ='direct_e_regulon_metadata_filtered',
    subset_eRegulons = ['TBX21','RORC'],
    subset_regions = hvr,
    subset_genes = hvg,
    add_differential_gene_expression = True,
    add_differential_region_accessibility = True,
    differential_variable = ['SampleID'])

G, pos, edge_tables, node_tables = create_nx_graph(nx_tables,
                   use_edge_tables = ['TF2R','R2G'],
                   color_edge_by = {'TF2R': {'variable' : 'TF', 'category_color' : {'TBX21':'Orange','RORC': 'Purple'}},
                                    'R2G': {'variable' : 'importance_x_rho', 'continuous_color' : 'Greys', 'v_min': -1, 'v_max': 1}},
                   transparency_edge_by =  {'R2G': {'variable' : 'importance_R2G', 'min_alpha': 0.1, 'v_min': 0}},
                   width_edge_by = {'R2G': {'variable' : 'importance_R2G', 'max_size' :  2, 'min_size' : 1}},
                   color_node_by = {'TF': {'variable': 'SampleID_Log2FC_Th1', 'continuous_color' : 'bwr','v_max':2, 'v_min':-2},
                                    'Gene': {'variable': 'SampleID_Log2FC_Th1', 'continuous_color' : 'bwr','v_max':5, 'v_min':-5},
                                    'Region': {'variable': 'SampleID_Log2FC_Th1', 'continuous_color' : 'PRGn', 'v_max':3,'v_min':-3}},
#                   transparency_node_by =  {'Region': {'variable' : 'CellType_Log2FC_Th1', 'min_alpha': 0.5},
#                                    'Gene': {'variable' : 'CellType_Log2FC_Th1', 'min_alpha': 0.5}},
                   size_node_by = {'TF': {'variable': 'fixed_size', 'fixed_size': 50},
                                    'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
                                    'Region': {'variable': 'fixed_size', 'fixed_size': 10}},
                   shape_node_by = {'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
                   label_size_by = {'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 20.0},
                                    'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 10.0},
                                    'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}},
                   layout='kamada_kawai_layout',
                   scale_position_by=250)

plt.figure(figsize=(10,10))
plot_networkx(G, pos)
plt.savefig('NCBR-160_scenicP_IL22IFNG_network_plot_redo.pdf')

export_to_cytoscape(G, pos, out_file = os.path.join('NCBR-160_scenicP_IL22IFNG_network_plot_redo.cyjs'))
