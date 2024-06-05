#!/usr/bin/env python3

# source .scenicPlus_env/bin/activate

import pandas as pd
import numpy as np
import pycisTopic
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.topic_binarization import *
from pycisTopic.topic_qc import *
from pycisTopic.diff_features import *
import pickle

# Load scATAC data
count_matrix2=pd.read_csv('peak_counts_Th1Polar.csv', sep=',',index_col=0)
cell_data =  pd.read_csv('ATAC_metadata_Th1Polar.txt', sep='\t',index_col=0)

# I chose "CellType" as the variable name to split the population on. 
#     The variable name needs to be constant between the scRNA and scATAC metadata files and have the same set of potential values.
# I needed to change the metadata information so that it would be identical between scATAC and scRNA
cell_data['CellType'] = np.where(cell_data['dataset'] == 'AB4710', 'Th1', 'Polar')

cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix2)
cistopic_obj.add_cell_data(cell_data)

# Load model data from pycistopic_model.py
models=pickle.load(open('Th1Polar_pycistopic_models_try1.pkl', 'rb'))
model = evaluate_models(models,
                        select_model=None,
                      return_model=True,
                      metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                      plot_metrics=False,
                      save = 'pycistopic_model_evaluation.pdf')

# Let the tool name the ideal model as seen in the pdf
cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj, open('pycistopic.pkl','wb'))

# As of June 2024, pycistopic has a new heatmap that may be informative, but it is not included in the docker container Skyler built.

# Find modules/topics of variable bins/peaks
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu', #ntop=3000,
                                         plot=True, num_columns=3, save= 'otsu.pdf')

binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li',
                                       plot=True, num_columns=3, nbins=100, save="cell_binarize.pdf")

topic_qc_metrics = compute_topic_metrics(cistopic_obj)

# Associate topics with the key metadata
topic_annot = topic_annotation(cistopic_obj, annot_var='CellType',
                               binarized_cell_topic=binarized_cell_topic, general_topic_thr = 0.2)

topic_qc_metrics = pd.concat([topic_annot[['CellType', 'Ratio_cells_in_topic',
                                           'Ratio_group_in_population','is_general']], topic_qc_metrics], axis=1)
topic_qc_metrics.to_csv('topic_annotation.csv')

# Find differentially accessible peaks
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None,
                                       selected_regions=None, scale_factor=10**6)

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='CellType',
                                  var_features=variable_regions, split_pattern = '-')

region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

pickle.dump(region_bin_topics_otsu, open('region_bin_topics_otsu.pkl', 'wb'))
pickle.dump(region_bin_topics_top3k, open('region_bin_topics_top3k.pkl', 'wb'))
pickle.dump(markers_dict, open('markers_dict.pkl', 'wb'))


