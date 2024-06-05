#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pycisTopic
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
import pickle

# Load scATAC data
count_matrix2=pd.read_csv('peak_counts.csv', sep=',',index_col=0)
cell_data =  pd.read_csv('ATAC_metadata.txt', sep='\t',index_col=0)

# I chose "CellType" as the variable name to split the population on. 
#     The variable name needs to be constant between the scRNA and scATAC metadata files and have the same set of potential values.
# I needed to change the metadata information so that it would be identical between scATAC and scRNA
cell_data['CellType'] = np.where(cell_data['dataset'] == 'AB4710', 'Th1', 'Polar')

cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix2)
cistopic_obj.add_cell_data(cell_data)

# Try various n_topics ranging between 2 and 50 until you get what you like
# n_cpu should be at least as many as the number of topics. Please do not
# set n_cpus to more than 8, it will cause ray to error out!
# Use "save_path" to save each model as separate pickles for later usage.
# NOTE: rerunning an individual model over and over will give the same result, but different combinations of models
#      may result in a different model chosen.
models=run_cgs_models(cistopic_obj,
                    n_topics=[5,10,15,20,25,30,40,50],
                    n_cpu=8,
                    n_iter=500, random_state=555,
                    alpha=50, alpha_by_topic=True,
                    eta=0.1, eta_by_topic=False,
                    save_path=None)

pickle.dump(models, open('pycistopic_models.pkl','wb'))

model = evaluate_models(models,
                        select_model=None,
                      return_model=True,
                      metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                      plot_metrics=False, # this errors out
                      save = 'pycistopic_model_evaluation.pdf')

