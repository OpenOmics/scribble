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
