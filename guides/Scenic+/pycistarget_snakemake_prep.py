#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
import pickle

# Load data from pycistopic_part2.py
region_bin_topics_otsu = pickle.load(open('region_bin_topics_otsu.pkl', 'rb'))
region_bin_topics_top3k = pickle.load(open('region_bin_topics_top3k.pkl', 'rb'))
markers_dict = pickle.load(open('markers_dict.pkl', 'rb'))

outFolder="bedfiles/"

# Extract peak coordinates for each topic
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')]
    #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    pr.PyRanges(region_names_to_coordinates(regions)).to_bed(path= outFolder + topic + "_topics_otsu_regions.bed")
    
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')]
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    pr.PyRanges(region_names_to_coordinates(regions)).to_bed(path=outFolder + topic + "_topics_top3k_regions.bed")

for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')]
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
    pr.PyRanges(region_names_to_coordinates(regions)).to_bed(path=outFolder + DAR + "_DAR_regions.bed")


# Warning messages will appear
