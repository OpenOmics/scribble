#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""fasta_filter.py: Filters a fasta file to keep sequences 
greater than or equal to a specified length.

Also, formats/folds the output so sequences do not exceed a
length of 80 characters. Any sequences exceeding this length
are wrapped to a new line. Output is directed to standard 
output and should be re-directed to a file.

Example usage:
$ ./fasta_filter.py /path/to/input.fa 1000 > seqs_greater_1000_length.fa
"""

from __future__ import print_function
import sys


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


def fold(sequence, max_length=80):
    """Folds or formats a long sequence so it does not exceed `max_length`.
    If a `sequence` does exceeed the limit, it will wrap sequenced onto a new line.
    @param sequence <str>:
        Coding DNA reference sequence to translate (transcript sequence or CDS sequence)
    @returns folded sequence <str>:
        Folded or formatted sequence that does not exceed `max_length`
    """
    # Clean sequence prior to conversion
    sequence = sequence.strip()

    # Format sequence to max length 
    chunks = []
    for i in range(0, len(sequence)-1, max_length):
        chunk = sequence[i:i+max_length]
        chunks.append(chunk)
    
    return "\n".join(chunks)


def main():
    """Main method, entry point that runs when program is directly invoked.
    """
    # Check for required positional arguments
    if len(sys.argv) == 1 or ('-h' in sys.argv or '--help' in sys.argv):
        print('USAGE: {} /path/to/input.fasta 1000 > output.fasta'.format(sys.argv[0]))
        sys.exit(1)

    # Input fasta file to filter
    input_fa = sys.argv[1]
    # Keep sequences greater than or
    # equal to this length
    try: filter_threshold = sys.argv[2]
    except IndexError: filter_threshold = 1000
    filter_threshold = int(filter_threshold)

    for sid, seq in fasta(input_fa):
        if len(seq) >= filter_threshold:
            print(
                "> {0}\n{1}".format(
                    sid, fold(seq, max_length=80)
                )
            )


if __name__ == '__main__':
    main()