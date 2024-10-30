# Oneliners

The onliners directory contains a collection of oneliners. These are helpful shell commands-- that may or may not make blow your mind. Did you write a oneliner that automates your last job? Did your perl oneliner regex talk back to you? Scribble it down here!

### For changing wall time for all jobs in the queue:

*What?* This oneliner will update the wall time of all of a user's jobs on Biowulf. In this example, we set the new wall time to 1 day.  

*Why?* This oneliner is useful for quickly increasing your walltimes. It can be useful for modifying your wall times when there is an approaching scheduled downtime for maintenance.

**Please note:** This oneliner will only work on Biowulf due to the use of the `newwall` command.


```bash
for ID in $(sjobs --tab --noheader | cut -f2); do newwall -j $ID -t 1-00:00:00; done
```


### Filter BED-like file for main reference chromosomes

*What?* This oneliner will remove non-mainline chromosomes from a BED-like file. Meaning, that it will only contain the following chromosomes after running: **chr1-22, chrM, chrX, and chrY**. It can be used with the following BED-like file formats: VCF, GTF, GFF, and BED. 

*Why?* This oneliner is useful for quickly removing scaffolds, decoys, and patches from a BED-like file. 

**Please note:** This example is more or less tailored for the human reference genome; however, it could easily be edited to be applied to other reference genomes with minimal editing-- such as a mouse reference genome or genomes from external providers that do not use the `chr` prefix for chromosomes. 

**TLDR:** The main takeaway from this example is that we can give grep a range or sequence of characters. The first pattern is looking for comments. The second pattern in the grep command is looking for lines that start with chr1-9 and chrM/Y/X. The third pattern is looking for lines that start with chr10-19, and the last pattern is looking for lines that start with chr20-22. Together this will capture everything outlined above. 

```bash
grep --color=never '^#\|^chr[1-9XYM]\b\|^chr1[0-9]\b\|^chr2[0-2]\b' /path/to/input.gtf > /path/to/output.gtf
```

### Find unique combinations of a set of files

*What?* This bash function will find the unique combinations of a set of files.

*Why?* This is useful for evaluating the results of different sets of tools. Let's say you would to understand the intersection of results from an ensembl of tools. This function will generate a space delimeted list of unique combinations of the input files.

```bash
function unique_combinations() {
    # Returns all unique possible combinations
    # of a collection of elements of a given length, 
    # to return all possible unique combinations, 
    # set the last argument equal to the length of 
    # the input collection, see example below: 
    # @Input $1 to (N-1): Set of items to permutate
    # @Input $N: Cardinality of the unique combination sets, can include a range or a set number
    # @Example: 
    # $ unique_combinations A B C 2-3
    # > A B
    # > A C
    # > B C
    # > A B C
    # Check if python is installed
    command -V python &> /dev/null || \
        { cat <<< 'Error: missing python, exiting now!' 1>&2; return 1; }
        
cat << EOF | python /dev/stdin "$@"
from __future__ import print_function
import sys

def unique_n_len_combinations(items, n):
    if n==0: yield []
    else:
      for i in range(len(items)):
            for cc in unique_n_len_combinations(items[i+1:], n-1):
                yield [str(items[i])] + cc

# Parse command line arguments
s = sys.argv[1:-1]
n = sys.argv[-1]
range_len = len([ _ for _ in n.split('-') if _ ])
# Check if range is specified
start = 0 if '-' not in n else int(n.split('-')[0])-1
stop = int(n) if '-' not in n \
else len(s) if range_len==1 else int(n.split('-')[1])
# Generate unique combinations
# differing set lengths
for i in range(start,stop):
    for comb in unique_n_len_combinations(s, i+1):
        print(' '.join(comb))
EOF
}
```