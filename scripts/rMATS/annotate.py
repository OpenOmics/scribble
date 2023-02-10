#!/usr/bin/env python
from __future__ import print_function
import sys

usage = \
"""
USAGE:
  $ ./annotate.py [-h] \\
      <{AS_Event}.MATS.JC.txt> \\
      <fromGTF.{AS_Event}.txt> \\
      <fromGTF.novelJunction.{AS_Event}.txt> \\
      <fromGTF.novelSpliceSite.{AS_Event}.txt>
    
SYNOPSIS:
  Annotates each recorded rMATS event as either a known,
  novelJunction, or novelSpliceSite alternative splicing
  event. This script add a new column called eventType as 
  the last column of the file. The output of this script 
  is redirected to standard output.

  Please note: the order of the positional arguments to 
  this script matters. For more information about each of
  rMATS output file(s), please see its documentation:
    https://github.com/Xinglab/rmats-turbo#output

EXAMPLE:
  $ ./annotate.py \\
      SE.MATS.JC.txt \\
      fromGTF.SE.txt \\
      fromGTF.novelJunction.SE.txt \\
      fromGTF.novelSpliceSite.SE.txt \\
    > SE.MATS.JC.annotated.txt
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


def add_event_type(file, event_dict, event_type):
    """Parses provided input file to add an alternative splicing
    event and its event type to a dictionary.
    @params file <str>:
        Path to file to parse
    @params event_dict <dict[id]=[event_type,..]>
        Dictionary to store each recorded alternative splicing event 
        with its type (i.e. known, novelJunction, novelSpliceSite), 
        where each key is the ID of the event and the value is a list
        containing strings of event types (1:M relationship).
    @params event_type <str>:
        The type of alternative splicing event: 
        known, novelJunction, novelSpliceSite.
    @returns updated_event_dict <dict[id]=[event_type,..]>:
        Updated events dictionary each new alternative splice sites type  
    """
    updated_event_dict = event_dict

    with open(file) as fh:
        # Skip over header 
        header = next(fh)
        for line in fh:
            linelist = line.lstrip().rstrip().split('\t')
            eid = linelist[0]
            if eid not in updated_event_dict:
                # The relationship between an event 
                # and its type many be one to many,
                # append to a list 
                updated_event_dict[eid] = []
            updated_event_dict[eid].append(event_type)
    
    return updated_event_dict


def main():
    # Check for help
    if '-h' in sys.argv or '--help' in sys.argv:
        err(usage)
        sys.exit()

    try:
        # Parse args
        # Check for usage, NOTE: the order 
        # of positional args matter
        rmats_results = sys.argv[1]
        all_splicing = sys.argv[2]
        novel_junctions = sys.argv[3]
        novel_splice_sites = sys.argv[4]
    except IndexError: 
        # Missing required args
        err('Error: missing one or more required argument(s)!')
        fatal(usage)

    # Dictionary to map event ID
    # to its type, i.e. whether it
    # is a novel or known alternative
    # splicing event
    index2type = {}
    # Annotate novel junctions
    index2type = add_event_type(
        file = novel_junctions,
        event_dict=index2type,
        event_type="novelJunction"
    )
    # Annotate novel splice sites
    index2type = add_event_type(
        file = novel_splice_sites,
        event_dict=index2type,
        event_type="novelSpliceSite"
    )
    # Add each alternative splice sites
    # type information, i.e. whether its
    # known or novel, etc.
    with open(rmats_results) as fh:
        header = next(fh)
        print("{0}\teventType".format(header.rstrip()))
        for line in fh:
            linelist = line.lstrip().rstrip().split('\t')
            eid = linelist[0]
            etype = 'known'
            if eid in index2type:
                # Add its event type information,
                # i.e. whether it is a novelJunction 
                # or a novelSpliceSite
                etype = ','.join(index2type[eid])
            linelist.append(etype)
            # Print added annotation information 
            # to standard output  
            print('\t'.join(linelist))


if __name__ == '__main__':
    main()