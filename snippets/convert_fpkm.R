library(edgeR)

# lines 6-9 are specific to the NCI Frederick sequencing facility outputs
# The idea is to make a matrix of the rsem data with the row.names being gene IDs
counts <- read.delim("RawCountFile_rsemgenes.txt", header = TRUE, stringsAsFactors = FALSE)
row.names(counts) <- counts$gene_id
counts <- counts[,-1]
counts <- counts[order(row.names(counts)),]

# geneinfo.bed comes from CCBR Pipeliner, chosen because the number of rows already match the number of rows in the rsem file
# the point is to make a data.frame with row.names in the equivalent format to the rsem matrix
# and the single column of data points being the gene lengths
gtf <- read.delim("geneinfo.bed", header = FALSE, stringsAsFactors = FALSE)
gtf$length <- (gtf$V3 - gtf$V2)
lengths <- data.frame(row.names = paste0(gtf$V5,"_",gtf$V7), length= gtf$length)

countslengths <- merge(counts, lengths, by=0, all=TRUE)
DGE <- DGEList(counts)
DGE <- calcNormFactors(DGE)
FPKM <- as.data.frame(rpkm(DGE, gene.length = countslengths$length, normalized.lib.sizes = TRUE))
write.table(FPKM,"FPKM.all_samples.txt",sep="\t",quote=F)
