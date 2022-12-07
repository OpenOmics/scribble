# Oneliners


For changing wall time for all jobs in the queue:  
``for ID in `sjobs --tab --noheader | cut -f 2`; do newwall -j $ID -t 1-00:00:00; done``
