# Name: merge_idxstat.R
# By: Tovah Markowitz
# Purpose: Have a bunch of idxstat output files and want to merge them into one large table?
# Use this snippet to load them all in, merge, and save in a txt file


inFiles <- list.files(pattern=".idxstat$")

sampleNames <- gsub(".idxstat","",inFiles)

inData <- read.delim(inFiles[1],header=F)
DataAll <- inData[,1:3]
names(DataAll) <- c("chr","length",sampleNames[1])

for (i in 2:length(inFiles)) {
  inData <- read.delim(inFiles[i],header=F)
  DataAll <- data.frame(DataAll,tmp=inData[,3])
  names(DataAll)[i+2] <- sampleNames[i]
}

write.table(DataAll, "idxstat_table.txt",quote=F,sep="\t",row.names=F)
