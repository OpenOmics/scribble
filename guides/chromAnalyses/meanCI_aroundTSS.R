####################
#
# Name: MeanCI_aroundTSS.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: December 6, 2023
#
####################
#
# Example code to:
# 1. Extract data from specific windows from multiple bigwigs.
#     In this case, we are extracting the information at TSS sites
#     and ordering all sites to have the first codon to the right of
#     of the TSS.
# 2. Calculate the mean across the these regions for the given sample.
# 3. Create a faceted plot of the results.
# 
# Note: this code also includes splitting the genes into groups based upon
#      external information about expected expression.

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(EnrichedHeatmap)

bws <- list.files(path="bigwig", pattern=".Q5DD.RPGC.inputnorm.bw", full.names=T)

# this is specifically for a file that is part of the cfChIP pipeline
# will also work with other txt files as long as the columns are labelled correctly
#    adapt as needed
# if running this on other samples with a bed file or gtf file annotation,
#    load with rtracklayer
PrebedFile <- "hg19.ensembl.prot_coding.with_annotations.txt"
inBed <- makeGRangesFromDataFrame( read.delim(PrebedFile), keep.extra.columns=T)

# define the gene starts and orientation in a way that EnrichedHeatmap can interpret
# The example here is for TSS sites, but it will work for any feature.
# Just adjust these three lines.
bed_starts <- inBed
start(bed_starts) <- ifelse( strand(inBed) == '+', start(inBed), end(inBed) )
end(bed_starts) <- start(bed_starts)

# this is if you want panels on the graph
# in this example code these strings were in the sample names
Groups <- c("E16","S16","E3","S3")

# this variable is also specific to the cfChIP txt file
# this is to split the genes into two groups on the same plot
Exp <- c("Housekeeping","NotExp") 

# this can take up to 10min per sample so there are multiple print statements
# to make sure the job is still progressing
for (j in 1:length(Groups)) {
  Group <- Groups[j]
  print(Group)
  bws2 <- grep(Group, bws, value=T)
  for ( i in 1:length(bws2) ) {
    print(i)
    signal <- import(bws2[i])
    for( k in 1:length(Exp)) {
      print(Exp[k])
      # grab only the relevant rows of the granges object for this mean line
      bed_starts2 <- bed_starts[which(data.frame(mcols(bed_starts)[Exp[k]]) == TRUE)]
      # adjust extend and w as needed where extend is the distance upstream and 
      # downstream of the granges window in order, and w is the window size to be
      # averaged. Use a higher window size to decrease the size of the final image
      # as files can become large but greater than 300 bins per inch on a graphic
      # tends to be useless
      tmp <- normalizeToMatrix(signal, bed_starts2,
                               extend=c(3000, 3000), w=25,
                               mean_mode="weighted", value_column="score")
      tmpMean <- colMeans(tmp, na.rm = T)
      # the next 7 lines are to allow all the data from the entire loop to be 
      # consolidated into a single data.frame. Adjust accordingly.
      # Note: Position is to help re-identify the TSS site for the downstream
      # figure and is different for a given set of extend/w parameters.
      tmpMean2 <- data.frame(Exp=Exp[k], Data=Group, SampleID=bws2[i],
                                     Position=seq(-119, 120), Mean=tmpMean)
      if( (i==1) & (j==1) & (k==1) ) {
          sepMeans <- tmpMean2
      } else {
          sepMeans <- rbind(sepMeans, tmpMean2)
      }
    }
  }
}


cols4 <- c("grey30","grey80",'#260F99','#8F7EE5','#6B990F','#C3E57E',
           '#990F0F','#FFB2B2','#99540F','#E5B17E')


ggplot(sepMeans, aes(x=Position, y=Mean, group=Exp, fill=Exp,colour=Exp))+
  stat_summary(geom="ribbon", fun.data=mean_cl_normal,
               fun.args=list(conf.int=0.75), alpha=0.2,color=NA) +
  stat_summary(geom="line", fun=mean,lwd=1) +
  facet_grid(Data ~ ., scales = "free_y") +
  labs(x = "", y = "Average RPGC") + theme_bw() +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(limits=c(-120,120),
                     breaks = c(-120, 0, 120),
                     labels = c("-3kb", "TSS","3kb")) +
  scale_colour_manual(values=cols4) +
  scale_fill_manual(values=cols4)
