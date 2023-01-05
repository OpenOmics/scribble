# Oneliners

The onliners directory contains a collection of oneliners. These are helpful shell commands-- that may or may not make blow your mind. Did you write a oneliner that automates your last job? Did your perl oneliner regex talk back to you? Scribble it down here!

### For changing wall time for all jobs in the queue:

*What?* This onliner will update the wall time of all of a user's jobs on Biowulf. In this example, we set the new wall time to 1 day.  

*Why?* This oneliner is useful for quikckly increasing your walltimes. It can be useful for modifying your wall times when there is an approaching scheduled downtime for maintenance.

**Please note:** this oneliner will only work on Biowulf due to the use of the `newwall` command.


```bash
for ID in $(sjobs --tab --noheader | cut -f2); do newwall -j $ID -t 1-00:00:00; done
```
