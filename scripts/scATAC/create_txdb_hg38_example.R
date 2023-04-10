library(GenomicFeatures)
library(AnnotationDbi)
library(data.table)

# Code used to create protein coding txdb for hg38
# Since nearly every line needs to be edited I left the code as is for an example.
# Notated lines are those to be edited.
# Please see https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.html for additional info.

chrom_file <- data.table(read.table('hg38.fa.sizes', col.names = c("chromosome", "size"))) # chromosome sizes file

seq.info <- Seqinfo(seqnames = chrom_file$chromosome, 
                    seqlengths = chrom_file$size, 
                    isCircular = c(as.vector(rep(FALSE, 24)), TRUE), # 24 chromosomes not circular.  MT is circular.
                    genome = "GRCh38") # Genome

txdb <- makeTxDbFromGFF("genes_protein_coding.gtf", # name of gtf, only protein coding were used for ATAC-seq
                        format = "gtf", 
                        dataSource = "ATAC-seq Pipeline GTF Gencode 28", # Notate gtf source
                        organism = "Homo sapiens", # organism
                        chrominfo = seq.info)

saveDb(txdb, "ATACseq_pipeline_protein-coding_hg38.txdb") # name of txdb to be created
