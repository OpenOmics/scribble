#!/usr/bin/env python3

"""
Name: CreatingIGVbatchFromBedFiles.py
Created by: Tovah Markowitz
Date: 9/14/18
Updated: 11/29/18

Purpose: To create an IGV batch script from a bed file or bed-like file
that will allow an semi-automated creation of snapshots across a genome. Just run
this script, open IGV, load the proper tracks and format them as you chose and then 
run the batch script to get a folder of pngs named by snapshop region.
"""
##########################################
# Modules
import optparse
import os

##########################################
# Functions

def read_bed(bedFile):
    # Purpose: read in bed file or any file with the first three columns being
    # in bed format and extract information
    f = open(bedFile, 'r')
    bed = f.readlines()
    f.close()
    bed2 = [ ]
    for i in range( len(bed) ):
        bed[i] = bed[i].strip().split('\t')
        if bed[i][0].startswith('chr'):
            bed2.append([bed[i][0]] + list(map(int, bed[i][1:3])) )
    return bed2

def create_IGV_batch(bedFile,bed,directory,offset=2000):
    # Purpose: to create an IGV batch script which says to load the bedfile,
    # where to save the pngs, and what regions to image.
    # It gives the option of defining how far on either side of the identified
    # regions to view.
    filename = bedFile.split(".")[0] + '_IGVbatch.txt'
    if not os.path.isdir(directory):
        os.mkdir(directory)
    f = open( filename, 'w')
    f.write('snapshotDirectory ' + directory + '\n')
    for row in bed:
        start = row[1] - int(offset)
        end = row[2] + int(offset)
        f.write('goto ' + row[0] + ':' + str(start) + '-' + str(end) + '\n' +
                        'snapshot' + '\n')
    f.close()

##########################################
# Main

def main():

    desc="""
    A script to create a simple batch script to make IGV images of multiple regions
    from a bed file (or bed-like file). To run the batch script, open IGV, load the
    genome-wide files of interest, and adjust tracks as needed. Then run the batch
    script (Tools -> Run Batch Script). Pngs will be created in chosen directory with 
    each image identified by genome position. Remember that like all other IGV snapshots
    the IGV window must remain unblocked/uncovered while the batch script is running.
    """
    
    parser = optparse.OptionParser(description=desc)
    
    parser.add_option('-b', dest= 'bedFile', default= '', help= "This is the name \
of the the input bed (or bed-like) file. Must be an absolute path.")
    parser.add_option('-d', dest= 'directory', default= '', help= "This is the \
name of the output directory where IGV should save the file. Must be an \
absolute path.")
    parser.add_option('-o', dest='offset', default= 2000, help="This is the \
distance to extend the visualizations on either side of a peak in base \
pairs. Default is 2000 bases.")

    (options,args) = parser.parse_args()
    
    bedFile = options.bedFile
    directory = options.directory
    offset = options.offset
    
    bed = read_bed(bedFile)
    create_IGV_batch(bedFile,bed,directory,offset)

if __name__ == '__main__':
    main()