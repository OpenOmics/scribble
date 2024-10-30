####################
#
# Name: igv_plot_from_r.R
# Created by: Subrata Paul, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: May 3, 2024
# Updated: May 3, 2024
# 
####################
#
# Purpose: To get igv like figure from r

##########################
# Inputs
# Path to tab seperated sample information file
# Required columns: file, name

samples<-'./AnalysisATAC/samples.txt'

# result file contain RNAseq DE and ATACseq DB analysis result in tsv format
# Required columns: pid - peak id in the format of chr*:from_to
# gene_name
# DB_logFC, DB_P  - differential binding logFC and P
# DE_logFC, DE_P - differential expression logFC and P

results<-'./temp/dream2.PvsHCinTNF.de.db.full.tsv'

# idoTrack file
# Sometimes getting information from different server on Biowulf is blocked
# The idoTrack for all chromosomes are saved locally
# It is in the .rds object contaiing each chromosome
idoTrack.file = '/data/NIAMS_IDSS/projects/NIAMS-42/Reports/Resources/all.idoTracks.rds'


#Use
# tracks = igvBWinR(samples, results, gene, idoTrack.file)
# Gviz::plotTracks(tracks, type = 'histogram', ylim = c(0,100),cex.title = 1, col.title = 'red')



igvBWinR<-function(samples, results, gene, idoTrack.file){
  if(is.character(samples)){
    if(file.exists(samples)){
      samples = read.table(samples, header = T)
    }else{
      stop('Sample file does not exist')
    }
  }
  
  if(!is.data.frame(samples)){
    stop('Samples must be a data.frame or file path')
  }
  
  if(is.character(results)){
    if(file.exists(results)){
      results = read.table(results, header = T)
    }else{
      stop('Result file does not exist')
    }
  }
  
  if(!is.data.frame(results)){
    stop('Result must be a data.frame or file path')
  }
  
  if(!all(c('pid', 'DB_logFC', 'DB_P', 'DE_logFC', 'DE_P')%in%names(results))){
    stop('Result file does not contains all required columns')
  }
  
  if(!all(c('file', 'name')%in%names(samples))){
    stop('Sample file does not contain required columns')
  }
  
  
  
  top.peaks = results%>%
    dplyr::filter(gene_name == gene)%>%
    dplyr::filter(DB_P<0.05)%>%
    dplyr::pull(pid)
  
  top.peaks.ranges = gsub('.*.:', '', top.peaks)%>%
    stringr::str_split(.,'_', simplify = T)%>%
    data.frame()%>%
    dplyr::mutate_all(as.numeric)
  
  all.idoTracks = readRDS(idoTrack.file)
  
  ensdb = filter(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, 
                 filter = ~ gene_name == gene)%>%
    getGeneRegionTrackForGviz()
  
  if(length(ensdb)==0){
    return(ggplot()+geom_text(aes(x = 0, y =0, label = paste0('Gene ', gene, ' was not found in EnsDb.Hsapiens.v86')))+theme_classic()+theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()))
    next
  }
  
  ensdb <- keepSeqlevels(ensdb, standardChromosomes(ensdb), pruning.mode = 'coarse')
  
  biomTrack = GeneRegionTrack(ensdb, geneSymbol = T, showId = T, transcriptAnnotation = 'transcript')
  idoTrack = all.idoTracks[[biomTrack@chromosome]]
  gtrack <- GenomeAxisTrack()
  
  # All tracks
  
  all_tracks <- list(idoTrack, gtrack)
  for(i in 1:nrow(samples)){
    all_tracks[[i+2]]<-DataTrack(samples$file[i],
                                 genome = biomTrack@genome,
                                 chromosome = biomTrack@chromosome,
                                 name = samples$name[i])
  }
  
  all_tracks[[length(all_tracks)+1]]<-biomTrack
  
  if(any(lapply(all_tracks, is.null)%>%unlist())){
    return(ggplot()+geom_text(aes(x = 0, y =0, label = paste0('One of the track was emptly for Gene ', gene)))+theme_classic()+theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()))
  }
  # Annotate peaks
  
  ht <- HighlightTrack(trackList = all_tracks,
                       start = top.peaks.ranges$X1, end = top.peaks.ranges$X2,
                       chromosome = 1, alpha = 0.2)  
  return(ht)
}

