## DEPRECIATED ##

#!/usr/bin/env python3

# source .scenicPlus_env/bin/activate

import pickle
import dill
from scenicplus.wrappers.run_scenicplus import run_scenicplus

scplus_obj = dill.load(open('scplus_obj.pkl', 'rb'))

# see scenic+_prep.py to chose the version
biomart_host = "http://sep2019.archive.ensembl.org/"
tf_file = "utoronto_human_tfs_v_1.01.txt"

# Make sure you have lots of free space in TMPDIR. It will use a lot for one of the steps and then blow everything away. I ended up needing 5T for a 5G .pkl file.
run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['CellType'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = tf_file,
        save_path = 'scenic',
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = False,
        export_to_loom_file = False,
        export_to_UCSC_file = False,
        path_bedToBigBed = 'bigbed',
        n_cpu = 1,
        _temp_dir = TMPDIR)
