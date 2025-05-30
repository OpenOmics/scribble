Purpose: To take an scATAC fragments file and convert it into pileup data in bigwig format to use outside of Signac

library(GenomicRanges)

sampleID <- "6_AB4834"
inFile <- paste0( sampleID, '/outs/fragments.tsv.gz')

outFile	<- paste0(sampleID, ".scATAC.bw")

inFrag <- read.delim(inFile, header=F)
colnames(inFrag) <- c("chr", "start", "end", "name")
fragDat <- makeGRangesFromDataFrame(inFrag)

covDat <- coverage(fragDat)
rtracklayer::export(covDat, outFile)
