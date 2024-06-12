##DEPRECIATED##

#!/usr/bin/env python3

# source .scenicPlus_env/bin/activate

import pandas as pd
import numpy as np
import scanpy as sc
import pickle
import dill
from scenicplus.scenicplus_class import create_SCENICPLUS_object

#################
# Step 1: Prep structure for scenic+

# Grab necessary scATAC information
menr = dill.load(open('motifs/menr.pkl', 'rb'))
cistopic_obj = pickle.load(open('pycistopic.pkl', 'rb'))

# Grab necessary scRNA data
# Confirm information is in the correct format with X holding the scaled.data and raw.X holding the data.
# https://mojaveazure.github.io/seurat-disk/reference/Convert.html
adata = sc.read_h5ad('scRNA2.h5ad')

# I needed to change the metadata information so that it would be identical between scATAC and scRNA
adata.obs['CellType'] = np.where(adata.obs['orig.ident'] == 'AB4720', 'Th1', 'Polar')

# change the nr_cells_per_metacells depending on the object size
# more cells for a larger object, less for a smaller one
# this example had ~20k scRNA and ~20k scATAC cells
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr,
    key_to_group_by = 'CellType',
    multi_ome_mode = False,
    meta_cell_split='.', # not '_' since that appears in cell type labels 
    nr_cells_per_metacells = 10
)

dill.dump(scplus_obj, open('scplus_obj.pkl', 'wb'))

#################
# Step 2: Find optimal host reference

ensembl_version_dict = {'105': 'http://www.ensembl.org',
                        '104': 'http://may2021.archive.ensembl.org/',
                        '103': 'http://feb2021.archive.ensembl.org/',
                        '102': 'http://nov2020.archive.ensembl.org/',
                        '101': 'http://aug2020.archive.ensembl.org/',
                        '100': 'http://apr2020.archive.ensembl.org/',
                        '99': 'http://jan2020.archive.ensembl.org/',
                        '98': 'http://sep2019.archive.ensembl.org/',
                        '97': 'http://jul2019.archive.ensembl.org/',
                        '96': 'http://apr2019.archive.ensembl.org/',
                        '95': 'http://jan2019.archive.ensembl.org/',
                        '94': 'http://oct2018.archive.ensembl.org/',
                        '93': 'http://jul2018.archive.ensembl.org/',
                        '92': 'http://apr2018.archive.ensembl.org/',
                        '91': 'http://dec2017.archive.ensembl.org/',
                        '90': 'http://aug2017.archive.ensembl.org/',
                        '89': 'http://may2017.archive.ensembl.org/',
                        '88': 'http://mar2017.archive.ensembl.org/',
                        '87': 'http://dec2016.archive.ensembl.org/',
                        '86': 'http://oct2016.archive.ensembl.org/',
                        '80': 'http://may2015.archive.ensembl.org/',
                        '77': 'http://oct2014.archive.ensembl.org/',
                        '75': 'http://feb2014.archive.ensembl.org/',
                        '54': 'http://may2009.archive.ensembl.org/'}

import pybiomart as pbm

def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
    print(f'host: {version}')
    try:
        n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
    except:
        print('Host not reachable')

v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")

