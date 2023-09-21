#!/usr/bin/env python3

import glob, sys, os

def main():
    files = glob.glob('*/outs/atac_fragments.tsv.gz')
    aggrFile = open('AggregatedDatasets.csv', 'w')

    #aggrFile = open('/Users/chenv3/Downloads/aggregate_samples.csv', 'w')

    aggrFile.write('library_id,atac_fragments,per_barcode_metrics,gex_molecule_info\n')
    files.sort()
    for i in files:
        sample = i.split('/')[0]
        if all([os.path.exists(f"{sample}/outs/per_barcode_metrics.csv"), os.path.exists(f"{sample}/outs/gex_molecule_info.h5")]):
            aggrFile.write(f"{sample},{os.path.abspath(i)},{os.path.join(os.path.abspath(sample), 'outs/per_barcode_metrics.csv')},{os.path.join(os.path.abspath(sample), 'outs/gex_molecule_info.h5')}\n")
    aggrFile.close()

if __name__ == "__main__":
    main()
