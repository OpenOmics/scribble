#!/usr/bin/env python3

import glob, sys, os

def main():
    fragment_files = glob.glob('*/outs/fragments.tsv.gz')
    cell_files = glob.glob('*/outs/singlecell.csv')

    samples = list(set([i.split('/')[0] for i in fragment_files]).intersection([i.split('/')[0] for i in cell_files]))
    samples.sort()
    aggrFile = open('AggregatedDatasets.csv', 'w')

    aggrFile.write('library_id,fragments,cells\n')
    for i in samples:
        aggrFile.write('%s,%s/outs/fragments.tsv.gz,%s/outs/singlecell.csv\n' % (i, os.path.abspath(i), os.path.abspath(i)))
    aggrFile.close()

if __name__ == "__main__":
    main()
