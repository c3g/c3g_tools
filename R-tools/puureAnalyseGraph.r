# echo "R --no-save --args 

# #########################################################################################################
# #########################################################################################################
# PARAM
# #########################################################################################################
# #########################################################################################################

args <- commandArgs(TRUE)
mainPath <- args[1]
mainPathSca <- args[2]
sample <- args[3]
type <- args[4]
fileRepeat <- args[5]
fileMap <- args[6]
fileAnno <- args[7]

algo <- ""
if (length(args)>7) {
  algo <- args[8]
}

#lib
library(grid)
library(ggplot2)

#set path 
folderScaffold <- paste(mainPath,"/scaffolds",algo,"/",sample,"/ray/ray21/",sep="")
folderOut <- paste(folderScaffold,"/insert",type,"/", sep="")
dir.create(file.path(folderOut), showWarnings = FALSE, recursive=TRUE)

# #########################################################################################################
# #########################################################################################################
# FONCTION
# #########################################################################################################
# #########################################################################################################


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


graphCov <- function(covIns, covSca, cluster, seqInsertion, sca, tabRepeatS, tabMapS, tabAnnoS, folderOut){

  listGGplot <- NULL

  if (is.null(cluster) | length(cluster[,1])==0){
    return()
  }
  
  minIns <- min(c(cluster$REF_START, cluster$REF_END))
  maxIns <- max(c(cluster$REF_START, cluster$REF_END))
  ggCovIns <- NULL
  if (length(covIns$ID_INS)>0){
    minIns <- min(covIns$POS)
    maxIns <- max(covIns$POS)
    ggCovIns <- ggplot(data=covIns, aes(x=POS, y=COV, group=TYPE, color=TYPE)) + geom_line(size=1.5)
    ggCovIns <- ggCovIns + theme(legend.background = element_rect(), legend.margin = unit(3.3, "cm")) + xlab(seqInsertion$CHR) +  ylab("Cov")
  }

  ggRepeat <- NULL
  if (length(tabRepeatS$CHR)>0){
    tabRepeatS$START[which(tabRepeatS$START < minIns)] <- minIns
    tabRepeatS$END[which(tabRepeatS$END > maxIns)] <- maxIns
    ggRepeat <- ggplot(data=tabRepeatS, aes(x=START, Y="R", color=TYPE)) + geom_segment(aes(y="R",yend="R",x=START,xend=END, color=TYPE)) + xlim(c(minIns, maxIns))
    ggRepeat <- ggRepeat + xlab("") +  ylab("Repeat") + theme(legend.background = element_rect(), legend.margin = unit(3.3, "cm"))
  }

  ggMap <- NULL
  if (length(tabMapS$CHR)>0){
    tabMapS$START[which(tabMapS$START < minIns)] <- minIns
    tabMapS$END[which(tabMapS$END > maxIns)] <- maxIns
    ggMap <- ggplot(data=tabMapS, aes(x=START, Y="M", color=TYPE)) + geom_segment(aes(y="M",yend="M",x=START,xend=END, color=TYPE)) + xlim(c(minIns, maxIns))
    ggMap <- ggMap + xlab("") +  ylab("Map") + theme(legend.background = element_rect(), legend.margin = unit(3.3, "cm"))
  }

  ggAnno <- NULL
  if (length(tabAnnoS$CHR)>0){
    tabAnnoS$START[which(tabAnnoS$START < minIns)] <- minIns
    tabAnnoS$END[which(tabAnnoS$END > maxIns)] <- maxIns
    ggAnno <- ggplot(data=tabAnnoS, aes(x=START, Y="A", color=TYPE)) + geom_segment(aes(y="A",yend="A",x=START,xend=END, color=TYPE)) + xlim(c(minIns, maxIns))
    ggAnno <- ggAnno + xlab("") +  ylab("Annot") + theme(legend.background = element_rect(), legend.margin = unit(3.3, "cm"))
  }

  ggCovSca <- NULL
  if (length(covSca$ID_SCA)>0){
    ggCovSca <- ggplot(data=covSca, aes(x=POS, y=COV, group=TYPE, color=TYPE)) + geom_line(size=1.5)
    ggCovSca <- ggCovSca + theme(legend.background = element_rect(), legend.margin = unit(3.3, "cm")) + xlab(paste("Scaffold ",sca$ID_SCA, " : ", sca$MATCH)) +  ylab("Cov")
  }

  tmp <- cluster$SCA_START[which(cluster$SCA_READ_FIRST==FALSE)]
  cluster$SCA_START[which(cluster$SCA_READ_FIRST==FALSE)] <- cluster$SCA_END[which(cluster$SCA_READ_FIRST==FALSE)]
  cluster$SCA_END[which(cluster$SCA_READ_FIRST==FALSE)] <- tmp

  tmp <- cluster$REF_START[which(cluster$SCA_READ_FIRST==FALSE)]
  cluster$REF_START[which(cluster$SCA_READ_FIRST==FALSE)] <- cluster$REF_END[which(cluster$SCA_READ_FIRST==FALSE)]
  cluster$REF_END[which(cluster$SCA_READ_FIRST==FALSE)] <- tmp
  tmp <- cluster$REF_START
  cluster$REF_START <- cluster$REF_END
  cluster$REF_END <- tmp

  ggClusterSca <- NULL
  ggClusterIns <- NULL
  if (length(cluster$REF_START)>0){
    
    ggClusterSca <- ggplot(data=cluster, aes(x=SCA_START, y=NUM_READ, group=TYPE, color=TYPE)) + xlim(c(min(covSca$POS), max(covSca$POS))) + ylim(c(0,max(cluster$NUM_READ))) + geom_segment(aes(xend=SCA_START+(SCA_END-SCA_START), yend=NUM_READ), size=1.5, alpha=.8, arrow=arrow(length=unit(0.1,"cm"))) + xlab("Scaffold seq") +  ylab("Num Sclip")
    ggClusterSca <- ggClusterSca + theme(legend.background = element_rect(), legend.margin = unit(3.3, "cm")) + xlab("") +  ylab("Num Read")

    ggClusterIns <- ggplot(data=cluster, aes(x=REF_START, y=NUM_READ, group=TYPE, color=TYPE)) + xlim(c(minIns, maxIns)) + ylim(c(0,max(cluster$NUM_READ))) + geom_segment(aes(xend=REF_START+(REF_END-REF_START), yend=NUM_READ), size=1.5, alpha=.8, arrow=arrow(length=unit(0.1,"cm"))) + xlab("Scaffold seq") +  ylab("Num Sclip")
    ggClusterIns <- ggClusterIns + theme(legend.background = element_rect(), legend.margin = unit(3.3, "cm")) + xlab("") +  ylab("Num Read")
  }

  if (!is.null(ggRepeat)){
    listGGplot <- c(listGGplot, list(ggRepeat))
  }
  if (!is.null(ggMap)){
    listGGplot <- c(listGGplot, list(ggMap))
  }
  if (!is.null(ggAnno)){
    listGGplot <- c(listGGplot, list(ggAnno))
  }
  listGGplot <- c(listGGplot,list(ggClusterIns, ggCovIns, ggCovSca, ggClusterSca))


  #if( !is.null(ggInsert[[1]]) | !is.null(ggInsert[[2]]) | !is.null(ggInsert[[3]]) ){
    #dir.create(file.path(paste(folderOut,"/",idSca,sep="")), showWarnings = FALSE, recursive=TRUE)
    h <- 12
    if(length(listGGplot)>4) {h <-14}
    if(length(listGGplot)>5) {h <-16}
    if(length(listGGplot)>6) {h <-18}
    pdf(file=paste(folderOut, "/", seqInsertion$ID_INS, ".pdf", sep=""), width=14, height=h)
    multiplot(plotlist=listGGplot, cols=1)
    dev.off()
  #}
}


# #########################################################################################################
# #########################################################################################################
# MAIN
# #########################################################################################################
# #########################################################################################################


file <- paste(folderOut,"scaffolds.filter.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
tabScaFilter <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character", "numeric", "numeric"))

file <- paste(folderOut,"insert.filter.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
allInsertFilter <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

file <- paste(folderOut,"cluster.filter.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
allClusterFilter <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "logical", "numeric", "character", "logical"))

file <- paste(folderOut,"covIns.filter.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
allCovInsFilter <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "numeric", "character"))

file <- paste(folderOut,"covSca.filter.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
allCovScaFilter <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "numeric", "character"))

file <- fileRepeat
if (!file.exists(file) | file.info(file)$size==0){
q()
}
tabRepeat <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=FALSE, colClasses=c("character", "numeric", "numeric", "character", "numeric", "character"))
tabRepeat <- cbind(tabRepeat, data.frame(TYPE=rep("Oth",length(tabRepeat[,1])), stringsAsFactors=FALSE))
select <- grep("Alu", tabRepeat[,4], fixed=TRUE)
tabRepeat$TYPE[select] <- "Alu"
select <- grep("L1", tabRepeat[,4], fixed=TRUE)
tabRepeat$TYPE[select] <- "L1"
select <- grep("LTR", tabRepeat[,4], fixed=TRUE)
tabRepeat$TYPE[select] <- "LTR"
select <- grep("MLT", tabRepeat[,4], fixed=TRUE)
tabRepeat$TYPE[select] <- "MLT"
select <- grep("HERV", tabRepeat[,4], fixed=TRUE)
tabRepeat$TYPE[select] <- "HERV"
select <- grep("MER", tabRepeat[,4], fixed=TRUE)
tabRepeat$TYPE[select] <- "MER"
colnames(tabRepeat) <- c("CHR","START", "END", "FAM", "LEN", "STRAND", "TYPE")
tabRepeat <- cbind(tabRepeat, data.frame(SOURCE=rep("Rmasker",length(tabRepeat[,1])), stringsAsFactors=FALSE))
tabRepeat <- tabRepeat[,c(1,2,3,7)]

file <- fileMap
if (!file.exists(file) | file.info(file)$size==0){
q()
}
tabMap <- read.table(file)
if (!is.numeric(tabMap[1,2])) {
tabMap <- read.table(file,header=TRUE)
}
colnames(tabMap) <- c("CHR","START", "END", "TYPE")
#tabMap$TYPE <- "Map"

file <- fileAnno
if (!file.exists(file) | file.info(file)$size==0){
q()
}
tabAnno <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=FALSE, colClasses=c("character", "numeric", "numeric", "character"))
colnames(tabAnno) <- c("CHR","START", "END", "TYPE")

for (i in 1:length(allInsertFilter$ID_INSERT)) {

  covIns <- allCovInsFilter[which(allCovInsFilter$ID_INS==allInsertFilter$ID_INSERT[i]),]
  covSca <- allCovScaFilter[which(allCovScaFilter$ID_SCA==allInsertFilter$ID_SCA[i]),]
  cluster <- allClusterFilter[which(allClusterFilter$ID_INSERT==allInsertFilter$ID_INSERT[i]),]

  seqInsertion <- allInsertFilter[i,]
  sca <- tabScaFilter[which(tabScaFilter$ID_SCA==allInsertFilter$ID_SCA[i]),]

  chr <- allInsertFilter$CHR[i]
  start <- min(covIns$POS)
  end <- max(covIns$POS)

  tabRepeatS <- tabRepeat[which( tabRepeat$CHR==chr& ((tabRepeat$START>=start & tabRepeat$START<=end) | (tabRepeat$END>=start & tabRepeat$END<=end)) ),]
  tabMapS <- tabMap[which( tabMap$CHR==chr & ((tabMap$START>=start & tabMap$START<=end) | (tabMap$END>=start & tabMap$END<=end)) ),]
  tabAnnoS <- tabAnno[which( tabAnno$CHR==chr & ((tabAnno$START>=start & tabAnno$START<=end) | (tabAnno$END>=start & tabAnno$END<=end)) ),]

  graphCov(covIns, covSca, cluster, seqInsertion, sca, tabRepeatS, tabMapS, tabAnnoS, folderOut)

}








