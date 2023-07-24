#!/usr/bin/env python3

# Python standard library
from __future__ import print_function, division
import os, sys

_help = """filter_results.py: Apply built-in filters rMATs results.
Usage:
    ./filter_results.py *.tsv 
"""


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


def index(linelist, columns):
    """Finds the index of each value listed in columns.
    @param linelist list[<str>]:
        Line split on its delimeter
    @params columns list[<str>]:
        Column names to index
    @return indices dict[<str>] = <int>:
        Dictionary containing the index of each col name
    """
    indices = {}
    missing = []
    for col in columns:
        try:
            i = linelist.index(col)
            indices[col] = i
        except ValueError as e:
            # Column is missing,
            # error out later to
            # see if anyother cols
            # are missing
            err("Error: Missing a required column... {}".format(col))
            missing.append(col)

    # Return errors if present
    if missing:
        fatal("Fatal: Missing the following required columns:\n\t-{}".format(
                missing 
            )
        )

    return indices


def avg(vals):
    """Find the average of a string of comma seperated values."""
    vlist = [int(v) for v in vals.split(',')]
    return sum(vlist) / len(vlist)


if __name__ == '__main__':
    # Input rMATs results
    try:
        input_files = sys.argv[1:]
    except IndexError:
        # No input file provided
        fatal(_help)

    # Create output directories
    flt1 = 'fdr_readcount'
    flt2 = 'fdr_readcount_incdiff'
    os.makedirs(flt1, exist_ok=True)
    os.makedirs(flt2, exist_ok=True)
    for file in input_files:
        # Add links to local images
        outfile = "{}.tsv".format(os.path.basename(file.rsplit('.',1)[0]))
        ofh1 = open(os.path.join(flt1, outfile), 'w')
        ofh2 = open(os.path.join(flt2, outfile), 'w')
        print('Filtering "{}" and writing to "{}"'.format(file, os.path.join(flt1, outfile)))
        print('Filtering "{}" and writing to "{}"'.format(file, os.path.join(flt2, outfile)))
        with open(file, 'r') as fh:
            # Index file header to extract
            # column information by name
            header = next(fh).rstrip().split('\t')
            new_header = header
            col2index = index(
                header, 
                ['FDR', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncLevelDifference']
            )
            # Add new columns for Filters
            ofh1.write("{}\n".format('\t'.join(new_header)))
            ofh2.write("{}\n".format('\t'.join(new_header)))
            for line in fh:
                # Apply the following filters, 
                # 1. FDR <= 0.05
                # 2. average(IJC_SAMPLE_1)>=10 || 
                #     average(SJC_SAMPLE_1)>=10 ||
                #     average(IJC_SAMPLE_2)>=10 ||
                #     average(SJC_SAMPLE_2')>-10
                # 3. abs(IncLevelDifference) >= 0.2
                linelist = line.rstrip().split('\t')
                fdr = float(linelist[col2index['FDR']])
                ijc1 = avg(linelist[col2index['IJC_SAMPLE_1']])
                sjc1 = avg(linelist[col2index['SJC_SAMPLE_1']])
                ijc2 = avg(linelist[col2index['IJC_SAMPLE_2']])
                sjc2 = avg(linelist[col2index['SJC_SAMPLE_2']])
                incdiff = float(linelist[col2index['IncLevelDifference']])
                if fdr <= 0.05 and (ijc1 >= 10 or sjc1 >= 10 or ijc2 >= 10 or sjc2 >= 10):
                    # print('\t'.join(linelist))
                    ofh1.write("{}\n".format('\t'.join(linelist)))
                    if abs(incdiff) >= 0.20:
                        ofh2.write("{}\n".format('\t'.join(linelist)))
        # Close output file handles
        ofh1.close()
