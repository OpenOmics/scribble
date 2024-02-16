#!/usr/bin/env python3

# source .scenicPlus_env/bin/activate

import pandas as pd
import numpy as np
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from scenicplus.wrappers.run_pycistarget import run_pycistarget
import pickle

# Load data from pycistopic_part2.py
region_bin_topics_otsu = pickle.load(open('region_bin_topics_otsu.pkl', 'rb'))
region_bin_topics_top3k = pickle.load(open('region_bin_topics_top3k.pkl', 'rb'))
markers_dict = pickle.load(open('markers_dict.pkl', 'rb'))

# Extract peak coordinates for each topic
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')]
    #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    pr.PyRanges(region_names_to_coordinates(regions)).to_bed(path=topic + "_topics_otsu_regions.bed")
    
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')]
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    pr.PyRanges(region_names_to_coordinates(regions)).to_bed(path=topic + "_topics_top3_regions.bed")

for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')]
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
    pr.PyRanges(region_names_to_coordinates(regions)).to_bed(path=DAR + "_DAR_regions.bed")

# Get required files from https://resources.aertslab.org/cistarget/
# These are the names of the ones needed for hg38
rankings_db = 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'
scores_db =  'hg38_screen_v10_clust.regions_vs_motifs.scores.feather'
motif_annotation = 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

# Don't forget to change the species name here
cistarget_dict = run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    save_path = 'motifs',
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    annotation = ['Direct_annot', 'Orthology_annot','Motif_similarity_annot',
                  'Motif_similarity_and_Orthology_annot'],
    n_cpu = 6,
    #_temp_dir = 'ray_spill',
    annotation_version = 'v10nr_clust',
    )

# All motif analysis results end up in save_path
# outputs are saved in menr.pkl
# html reports of all motif analyses will also be located there
