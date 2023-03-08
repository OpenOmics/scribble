#Single-cell ATAC R functions
#To use within Seurat/Signac framework

find_da_peaks <- function(sample_1, sample_2, seurat_obj, grouping, min_pct=0.1) {
  # find markers for given contrast
  # library(signac)
  # library(seurat)
  da_peaks <- FindMarkers(object = seurat_obj,
                                            ident.1 = sample_1,
                                            ident.2 = sample_2,
                                            test.use = 'LR',
                                            latent.vars = 'nCount_ATAC', # AKA nCount_peaks
                                            group.by = grouping,
                                            min.pct = min_pct) # Some cases will require 0.05
  file_name <- paste0(gsub(' ', '_', sample_1), "_vs_", gsub(' ', '_', sample_2), "_da_peaks.csv")
  write.csv(cbind(Peak = rownames(da_peaks), da_peaks), file_name, quote = FALSE, row.names = FALSE)

  return(da_peaks)
}

##################
# for annotations with ChIPseeker

# written for hg38, but easily adaptable

txdb <- loadDb("ATACseq_pipeline_protein-coding_hg38.txdb")

annotate_peaks <- function(seurat_object, da_peaks) {
  all_peaks_granges <- granges(seurat_object)
  da_peaks_granges <- granges(seurat_object[rownames(da_peaks),])

  peakAnnoList <- suppressWarnings(lapply(list(all_peaks_granges, da_peaks_granges), 
                       annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 2000), verbose=FALSE, annoDb="org.Hs.eg.db"))
  names(peakAnnoList) <- c('All Peaks', 'DA Peaks')
  return(peakAnnoList)
}

clean_DA_annotated <- function(peakAnnoList, da_peaks, sample_1, sample_2) {
  # Grab DA annotated peaks
  chipseeker_da_ann <- as.data.frame(peakAnnoList$`DA Peaks`)
  chipseeker_da_ann$annotation <- gsub(',', ';', chipseeker_da_ann$annotation)
  chipseeker_da_ann$geneChr <- paste0("chr", chipseeker_da_ann$geneChr)
  chipseeker_da_ann <- cbind(Peak = paste0(chipseeker_da_ann$seqnames, "-", chipseeker_da_ann$start,
                             "-", chipseeker_da_ann$end), chipseeker_da_ann)
  # write DA annotated peaks with p-values
   da_full <- merge(cbind(Peak = rownames(da_peaks), da_peaks), chipseeker_da_ann,
                    by.x = "Peak", by.y = "Peak")
   da_full <- da_full[order(match(da_full$Peak, rownames(da_peaks))),]
   file_name <- paste0(gsub(' ', '_', sample_1), "_vs_", gsub(' ', '_', sample_2),
                       "_da_peaks_ChIPseeker_ann_merged.csv")
   write.csv(da_full, file_name, quote = FALSE, row.names = FALSE)
   return(da_full)
}

##################
# motif analysis

# written for hg38 and JASPAR2022, but easily adaptable
# requires: TFBSTools

motif_prep <- function(seurat_object, txdb) {
  main.chroms <- standardChromosomes(txdb)
  keep.peaks <- as.logical(seqnames(granges(combined)) %in% main.chroms)
  seurat_object <- seurat_object[keep.peaks, ]

  pfm <- getMatrixSet( x = JASPAR2022,
    opts = list(collection = "CORE", species = 9606, all_versions = FALSE) )

  seurat_object <- AddMotifs(object = seurat_object,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm )
  return(seurat_object)
}

da_motifs <- function(seurat_object, da_peaks, sample_1, sample_2) {
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
  enriched.motifs <- FindMotifs(object = seurat_object,
                                features = top.da.peak)
   file_name <- paste0(gsub(' ', '_', sample_1), "_vs_", gsub(' ', '_', sample_2),
		       "_da_peaks_enriched_motifs.csv")
   write.csv(enriched.motifs, file_name, quote = FALSE, row.names = FALSE)
   return(enriched.motifs)
}

##################
# for interactive volcano plot
# requires plotly
hline <- function(y = 0, color = "black") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, width = 1, dash = "dot")
  )
}

vline <- function(x = 0, color = "black") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, dash="dot", width = 1)
  )
}

create_volcano_plotly <- function(sample_1, sample_2, da_peaks) {
    # table adjustments
  tmp <- da_peaks[,c(2,5)]
  tmp[which(tmp[,2] == 0),2] <- 1e-320
  
  # For color coding
  tmp$diffexpressed <- "No DE"
  tmp$diffexpressed[tmp$avg_log2FC > 0.5 & tmp$p_val_adj < 0.05] <- "Up"
  tmp$diffexpressed[tmp$avg_log2FC < -0.5 & tmp$p_val_adj < 0.05] <- "Down"
  
  volcano_colors <- vector()
  if (any(tmp$diffexpressed == "Down")) {
    volcano_colors <- c(volcano_colors, "red")
  }
  if (any(tmp$diffexpressed == "No DE")) {
    volcano_colors <- c(volcano_colors, "black")
  }
  if (any(tmp$diffexpressed == "Up")) {
    volcano_colors <- c(volcano_colors, "blue")
  }
  
  # labels for volcano plot
  tmp$delabel <- rownames(tmp)
  tmp$diffexpressed <- factor(tmp$diffexpressed, levels = c("Down", "No DE", "Up"))
  
  tmp$log_p <- -log10(tmp$p_val_adj)
  tmp <- tibble(tmp)
  
  volcano <- plot_ly(data = tmp, x = ~avg_log2FC, y = ~log_p, text = ~delabel, 
                     type="scatter", mode="markers", marker=list(size=4), color=~diffexpressed, 
                     colors=volcano_colors) %>%
    layout(shapes = list(vline(-0.5), vline(0.5), hline(-log10(0.05)))) %>%
    layout(title = paste0("DA Peaks: ", sample_1, " versus ", sample_2), 
           xaxis=list(title="log2 fold change"), yaxis=list(title="-log10 adjusted p-value"))
  volcano
}

EnhancedVolcano(da_full, lab = da_full$SYMBOL,
                x = 'avg_log2FC',
                y = 'p_val_adj', ylab = bquote(~Log[10]~ 'adj P'),
                drawConnectors = T, labSize= 3,
                title = paste0(group1, " vs ", group2))
