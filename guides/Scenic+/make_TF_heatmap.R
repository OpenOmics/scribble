# This is example code to grab the necessary files from the scenic plus outputs and create the TF heatmap shown in their papers/reports using R.

library(ggplot2)

outputFolder1 <- "scplus_CD4"
outputFolder2 <- "scplus_CD4/outs/"

eRegulon <- read.delim(paste0(outputFolder,"eRegulon_direct.tsv")) 
eReg <- unique(eRegulon[,c("TF", "Gene_signature_name", "Region_signature_name",
						   "eRegulon_name", "regulation")])

scaleDat <- function(x) {
  tmp <- as.matrix(x[2:nrow(x),2:ncol(x)])
  tmp <- apply(tmp,2,as.numeric)
  tmp <- apply(tmp,2,scale)
}


files <- list.files(path=outputFolder1, recursive=T, full.names=T, pattern="csv")

TF2G <- read.csv(files[1],header=F)
R2G <- read.csv(files[2],header=F)

# This was written for my analysis where there were only two conditions

TF2Ga <- data.frame(condition=c(rep(TF2G[2,1],ncol(TF2G)-1), rep(TF2G[3,1],ncol(TF2G)-1)),
                    Gene_signature_name=rep(unlist(TF2G[1,2:ncol(TF2G)]),2),
		    TF2G=c(unlist(TF2G[2,2:ncol(TF2G)]), unlist(TF2G[3,2:ncol(TF2G)])))

R2Ga <- data.frame(condition=c(rep(R2G[2,1],ncol(R2G)-1), rep(R2G[3,1],ncol(R2G)-1)),
                    Region_signature_name=rep(unlist(R2G[1,2:ncol(R2G)]),2),
                    R2G=c(unlist(R2G[2,2:ncol(R2G)]), unlist(R2G[3,2:ncol(R2G)])))

# scaled versions
TF2Gb <- scaleDat(TF2G)
R2Gb <- scaleDat(R2G)
TF2Ga <- data.frame(condition=c(rep(TF2G[2,1],ncol(TF2G)-1), rep(TF2G[3,1],ncol(TF2G)-1)),
                    Gene_signature_name=rep(unlist(TF2G[1,2:ncol(TF2G)]),2),
                    TF2G=c(TF2Gb[1,],TF2Gb[2,]))

R2Ga <- data.frame(condition=c(rep(R2G[2,1],ncol(R2G)-1), rep(R2G[3,1],ncol(R2G)-1)),
                    Region_signature_name=rep(unlist(R2G[1,2:ncol(R2G)]),2),
                    R2G=c(R2Gb[1,],R2Gb[2,]))



tmp <- merge(eReg,TF2Ga)
allData <- merge(tmp,R2Ga)
allData$regulation <- ifelse(allData$regulation == 1, "Activator", "Repressor")
allData$TF2G <- as.numeric(allData$TF2G)
allData$R2G <-	as.numeric(allData$R2G)
allData2 <- allData[which(endsWith(allData$eRegulon_name,"+")),]

p <- ggplot(allData2, aes(x=condition, y=eRegulon_name))
p <- p + facet_grid(regulation ~ ., scales="free", space = "free_y") +
    geom_tile(mapping=aes(fill=TF2G)) +
    scale_fill_distiller(type = 'div', palette = 'RdYlBu') +
    geom_point(mapping=aes(size=R2G), color = "black") +
    theme(axis_title_x = element_blank(), axis_title_y = element_blank()) +
    theme_bw()

