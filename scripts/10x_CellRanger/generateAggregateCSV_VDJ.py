#!/usr/bin/env python3

import glob, sys, os

def main():
    files = glob.glob('*/outs/vdj_contig_info.pb')
    aggrFile = open('AggregatedDatasets.csv', 'w')

    #aggrFile = open('/Users/chenv3/Downloads/aggregate_samples.csv', 'w')

    #aggrFile.write('sample_id,vdj_contig_info\n')
    aggrFile.write('library_id,vdj_contig_info,donor,origin\n')
    files.sort()
    for i in files:
        aggrFile.write('%s,%s\n' % (i.split('/')[0], os.path.abspath(i)))
    aggrFile.close()

if __name__ == "__main__":
    main()
