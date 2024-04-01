#!/usr/bin/env python3
"""
Parsing script for miRNA to RNA pairwise alignment output from
the miRanda v1.0b application output.

The miRanda output this works with is the std output captured in
the image: /data/OpenOmics/SIFs/miranda_1.0b.sif.

When run like: miranda <mature_miRNA> <transcriptome>

Uses pandas, stdlib [os, pathlib, regex, argparse]
"""
import re
import os
from pathlib import Path
import pandas as pd
import argparse


regex = r"^>{1,2}([a-zA-Z0-9-]+)\s+([A-Z_.0-9]+)\s+([0-9.]+)\s+([0-9-.]+)\s+([0-9-.]+)\s+([0-9-.]+)\s+([0-9-.]+)\s+([0-9-.]+)\s+([0-9-.]+)\s+([0-9-.]+)\s+([0-9]+)"


def is_valid_file(parser, arg):
    _p = Path(arg)
    if not _p.exists():
        parser.error("The file %s does not exist!" % arg)
    else:
        return str(_p.resolve())  # return an open file handle


def parse_text_block(text_block):
    matches = re.match(regex, text_block, re.MULTILINE)

    struct_data = None
    if matches:
        grps = matches.groups()
        assert len(grps) == 11, 'Invalid data lengths'
        struct_data = {
            'Seq1': grps[0],
            'Seq2': grps[1],
            'Tot Score': grps[2],
            'Tot Energy': grps[3],
            'Max Score': grps[4],
            'Max Energy': grps[5],
            'Strand': grps[6],
            'Len1': grps[7],
            'Len2': grps[8],
            'Positions': grps[9] + ':' + grps[10]
        }

    return struct_data


def main(args):
    parsed_data = None
    with open(args.input) as fi:
        lines = fi.readlines()
        # Parse text block
        parsed_data = list(filter(None, map(parse_text_block, lines)))

    if parsed_data:
        df = pd.DataFrame(parsed_data)
        df = df[df['Max Score'].str.replace('.', '').astype(float) > 0.0]
        df = df.sort_values(by=['Seq1', 'Seq2', 'Tot Score'], ascending=[False, False, True])
        df.to_csv(args.output, sep='\t', index=False)
    else:
        raise ValueError('Empty `parsed_data`: ' + str(parsed_data))
    

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=lambda x: is_valid_file(p, x), help='Input miRanda results in plain text')
    p.add_argument("output", type=lambda x: str(Path(x).resolve()), help='Output file name for parsed information')
    main(p.parse_args())
