#!/usr/bin/env python3

"""md5.py: calculates md5s of multiple files in parallel.

The md5 calculation is memory safe. It reads in a file
in blocks of 64 KiB.

USAGE:
  python3 md5.py file1.txt file2.csv --threads 8

Use --help to see all available options.
"""

# Python standard library
from __future__ import print_function
import argparse
import hashlib
import sys

# 3rd party imports
import ray


@ray.remote
def md5sum(filename, first_block_only = False, blocksize = 65536):
    """Gets md5checksum of a file in memory-safe manner.
    The file is read in blocks/chunks defined by the blocksize parameter. This is 
    a safer option to reading the entire file into memory if the file is very large.
    @param filename <str>:
        Input file on local filesystem to find md5 checksum
    @param first_block_only <bool>:
        Calculate md5 checksum of the first block/chunk only
    @param blocksize <int>:
        Blocksize of reading N chunks of data to reduce memory profile
    @return (filename, hasher.hexdigest()) <tuple>:
        filename and MD5 checksum of the file's contents
    """
    hasher = hashlib.md5()
    with open(filename, 'rb') as fh:
        buf = fh.read(blocksize)
        if first_block_only:
            # Calculate MD5 of first block or chunck of file.
            # This is a useful heuristic for when potentially 
            # calculating an MD5 checksum of thousand or 
            # millions of file.
            hasher.update(buf)
            return (filename, hasher.hexdigest())
        while len(buf) > 0:
            # Calculate MD5 checksum of entire file
            hasher.update(buf)
            buf = fh.read(blocksize)

    return (filename, hasher.hexdigest())


def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Calculates md5 checksums of multiple files in parallel.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'files',
        metavar='FILE',
        nargs='+',
        help='One or more files to checksum.',
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        help='Number of concurrent MD5 processes to spawn using ray.',
    )
    parser.add_argument(
        '-f', '--first-block-only',
        action='store_true',
        help='Only checksum the first block of each file (a fast heuristic '
             'when processing very large numbers of files).',
    )
    parser.add_argument(
        '-b', '--blocksize',
        type=int,
        default=65536,
        help='Number of bytes to read per chunk (reduces memory profile).',
    )
    parser.add_argument(
        '-d', '--delimiter',
        choices=('space', 'comma', 'tab'),
        default='space',
        help="Delimiter separating the fields in the output. 'space' (two "
             "spaces) matches the output format of the md5sum command.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # Map the delimiter choice to its actual character(s). The default,
    # two spaces, matches the output format of the md5sum command.
    delimiters = {'space': '  ', 'comma': ',', 'tab': '\t'}
    delimiter = delimiters[args.delimiter]

    # Initialize a ray cluster
    # with X remote workers
    ray.init(num_cpus = args.threads)

    # Run md5sum function in
    # parallel with the .remote()
    # method. This methods yields
    # a future or ObjectRef which
    # can be fetched later with
    # ray.get(). As files are
    # processed, remote workers
    # freed from a job queue.
    workers = [
        md5sum.remote(
            f,
            first_block_only=args.first_block_only,
            blocksize=args.blocksize,
        )
        for f in args.files
    ]

    # The ray.get() method is
    # blocking waits until all
    # tasks have completed. The
    # order of the inputs and
    # the function results are
    # preserved. Print out results
    # to standard output.
    for worker in workers:
        try:
            result = ray.get(worker)
            file, md5 = result
        except Exception as e:
            print('An error occured:\n{0}'.format(e), file=sys.stderr)
            continue

        print("{0}{1}{2}".format(md5, delimiter, file))


if __name__ == '__main__':
    main()
