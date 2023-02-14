#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import print_function
import sys

usage = \
"""
USAGE:
  $ ./rmats2VEPinput.py [-h] \\
      <genome.fa> \\
      <{AS_Event}.MATS.JC.filtered.txt> \\
      <{AS_Event}> \\
      <exonStart_0base_Column_Index> \\
      <exonEnd_Column_Index> \\
      <rMATS_{AS_Event}>
    
SYNOPSIS:
  Predicts the effect of a reported rMATS event 
  using VEP. Each AS event type is predicted as 
  a deletion. This script current only accepts 
  the following AS event types-- SE, A5SS, A3SS,
  and MXE-- where each AS event is treated as a 
  deletion.

  This script will coverts each recorded rMATS event
  into a format compatiable with VEP. Each AS event 
  is converted into the default VEP input (i.e. 
  Ensembl default). More information about VEP's 
  expected input format can be found here: 
  https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html

  PLEASE NOTE: the column indices for the 4th and 5th 
  positional arguments are zero-based. Please see the 
  examples below for running this script with the output 
  from different event types. There are example for SE
  and MXE events detected by rMATS.

EXAMPLE:
  # Create input for VEP using skipped
  # exon (SE) events detected by rMATS,
  # focusing on the effect of deleting
  # the reported skipped exon
  $ ./rmats2VEPinput.py \\
        GRCh38.primary_assembly.genome.fa \\
        SE.CNOT3.MATS.JC.filtered.txt \\
        SE \\
        5 \\
        6 \\
        rMATS_SE

  # Create input for VEP using mutually 
  # exclusive exon (MXE) events detected 
  # by rMATS, focusing on the effect of 
  # deleting the 1stExon in the pair.
  $ ./rmats2VEPinput.py \\
        GRCh38.primary_assembly.genome.fa \\
        MXE.CNOT3.MATS.JC.filtered.txt \\
        MXE \\
        5 \\
        6 \\
        rMATS_MXE_1stExon

  # Create input for VEP using mutually 
  # exclusive exon (MXE) events detected 
  # by rMATS, focusing on the effect of 
  # deleting the 2ndExon in the pair.
  $ ./rmats2VEPinput.py \\
        GRCh38.primary_assembly.genome.fa \\
        MXE.CNOT3.MATS.JC.filtered.txt \\
        MXE \\
        7 \\
        8 \\
        rMATS_MXE_2ndExon
"""


 # Helper Functions 
def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions 
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)


def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def fasta(filename):
    """
    Reads in a FASTA file and yields each of its entries.
    The generator yields each sequence identifier and its 
    corresponding sequence to ensure a low memory profile. 
    If a sequence occurs over multiple lines, the yielded 
    sequence is concatenated.
     @param filename <str>:
        Path of FASTA file to read and parse
    @yield chrom, sequence <str>, <str>:
        Yields each seq id and seq in the FASTA file
    """
    with open(filename, 'r') as file:
        sequence, chrom = '', ''
        for line in file:
            line = line.strip()
            if line.startswith('>') and sequence:
                # base case for additional entries
                yield chrom, sequence
                chrom = line[1:] # remove the > symbol
                sequence = ''
            elif line.startswith('>'):
                # base case for first entry in fasta file
                chrom = line[1:] # remove the > symbol
            else:
                # concatenate multi-line sequences
                sequence += line
        else:
            yield chrom, sequence


def revcomp(seq):
    """ Finds the reverse complement of a given sequence in upper case, 
    i.e a sequence on the '-' strand.
     @param seq <str>: Sequence to find reverse complement
     @returns rc <str>: Reverse complement of `seq` in upper case
    """
    complements = {
        "A": "T", "a": "T",
        "T": "A", "t": "A",
        "C": "G", "c": "G",
        "G": "C", "g": "C",
        "N": "N", "n": "N"
    }
    
    # Reverse the sequence
    seq = seq[::-1]
    rc = ''
    for bp in seq:
        try:
            comp = complements[bp]
        except KeyError as e:
            sys.exit(
                '{}\n\nFatal: encountered an non-standard nucleotide: {} @ {}'.format(
                    e,
                    bp,
                    seq
                )
            )
        # Build reverse complement string
        rc += comp
    
    return rc


def main():
    """
    Pseudo main method that runs when program is directly invoked.
    """
    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        err(usage)
        sys.exit()

    # Parse args and check for usage
    try:
        file = sys.argv[1]
        bed = sys.argv[2]
        as_type = sys.argv[3].upper()
        start_index = int(sys.argv[4])
        stop_index = int(sys.argv[5])
        vep_id_prefix = sys.argv[6]
    except IndexError:
        # Missing required args
        err('Error: missing one or more required positional argument(s)!')
        fatal(usage)

    # Create index of chrom seqs
    chrom2seq = {}
    for chrom, seq in fasta(file):
        chrom = chrom.split()[0]
        chrom2seq[chrom] = seq

    # Extract the sequence at a specific 
    # chrom, start, stop position 
    with open(bed, 'r') as fh:
        _header = next(fh)
        for line in fh:
            linelist = line.lstrip().rstrip().split('\t')
            # rMATS exon start is already 0-based
            # so we do not need to offset the start,
            # since the start is offset already the 
            # stop position is alerady inclusive, I 
            # confirmed this will a known AS event in
            # the results and compared it to the 
            # annotated exon in the GTF file
            eid    = str(linelist[0])             # rMATS event ID 
            geneID = str(linelist[1])             # GeneID
            gene   = str(linelist[2])             # geneSymbol
            start  = int(linelist[start_index])   # exonStart_0base
            stop   = int(linelist[stop_index])    # exonEnd
            chrom  = linelist[3].strip()          # chr
            strand = linelist[4].strip()          # strand
            seq    = chrom2seq[chrom][start:stop] # sequence of feature
            # Output the results in VEP default output format:
            # https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html
            print(
                "{0}\t{1}\t{2}\t{3}/-\t{4}\t{5}".format(
                    chrom.replace('chr', ''), 
                    start+1, 
                    stop, 
                    seq,
                    strand,
                    # VEP identifer contains,
                    # the event type (SE vs. MXE),
                    # the event ID from rMATS,
                    # the ensembl ID, and
                    # the gene symbol 
                    '{0}_{1}_{2}_{3}'.format(vep_id_prefix, eid, geneID, gene)
                )
            )

if __name__ == '__main__':
    main()
