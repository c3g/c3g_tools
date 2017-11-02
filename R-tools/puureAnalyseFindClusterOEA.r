# echo "R --no-save --args 
#‘TMPDIR’, ‘TMP’ and ‘TEMP’

# #########################################################################################################
# #########################################################################################################
# PARAM
# #########################################################################################################
# #########################################################################################################

args <- commandArgs(TRUE)
mainPath <- args[1]
kmer <- args[2]
sample <- args[3]
type <- args[4]
minMapq <- as.numeric(args[5])
insertSize <- as.numeric(args[6])

#lib
library(Rsamtools)
library(grid)
library(igraph) # install.packages("igraph", lib="~/.R/library")
library(GenomicAlignments)

#set path 
folderSC <- paste(mainPath,"/sclip/",sample,"/",sep="")
folderExtract <- paste(mainPath,"/extract/",sample,"/",sep="")
folderScaffold <- paste(mainPath,"/scaffolds/", sample,"/ray/ray",kmer,"/",sep="")
folderOut <- paste(folderScaffold,"/insert",type,"/", sep="")
dir.create(file.path(folderOut), showWarnings = FALSE, recursive=TRUE)


# #########################################################################################################
# #########################################################################################################
# GLOBAL
# #########################################################################################################
# #########################################################################################################
flagFirstRead<- scanBamFlag(isFirstMateRead=TRUE)
flagSecondRead<- scanBamFlag(isFirstMateRead=FALSE)

paramFirst <- ScanBamParam(flag=flagFirstRead)
paramSecond <- ScanBamParam(flag=flagSecondRead)

Gbam <- list(readGAlignments(paste(folderExtract,sample,".OEAMAP.bam",sep=""),param=paramFirst,use.names=TRUE), readGAlignments(paste(folderExtract,sample,".OEAMAP.bam",sep=""),param=paramSecond,use.names=TRUE))

# #########################################################################################################
# #########################################################################################################
# FONCTION
# #########################################################################################################
# #########################################################################################################

creatInsertCluster <- function(tab, idSca, type, first) {
  # crée les cluster ou les read map se chevauchent sur la ref
  gr <- GRanges(seqnames=tab$CHR, ranges = IRanges(start=tab$REF_START, end=tab$REF_END, names=tab$NAME))
  g <- graph.edgelist(as.matrix(findOverlaps(gr)), directed=TRUE)
  tabClusterRef <- data.frame(cluster=as.numeric(clusters(g)$membership), range=as.numeric(V(g)), key=paste(as.numeric(clusters(g)$membership),as.numeric(V(g)),sep="-"), stringsAsFactors=FALSE)

  # crée les cluster ou les read unmap se chevauchent sur le scaffolds
  gr <- GRanges(seqnames=rep(idSca,length(tab$SCA_START)), ranges = IRanges(start=tab$SCA_START, end=tab$SCA_END))
  g <- graph.edgelist(as.matrix(findOverlaps(gr)), directed=TRUE)
  tabClusterSca <- data.frame(cluster=as.numeric(clusters(g)$membership), range=as.numeric(V(g)), key=paste(as.numeric(clusters(g)$membership),as.numeric(V(g)),sep="-"), stringsAsFactors=FALSE)

  rm(gr)
  rm(g)
  #gc()

  # crée les cluster ou les reads se chevauchent sur le scaffolds et sur la ref
  cluster <- NULL
  for (i in unique(tabClusterRef$cluster)){
    ligneR <- tabClusterRef[which(tabClusterRef$cluster==i),]
    ligneR <- ligneR$range
    while(length(ligneR)!=0){
      clus <- tabClusterSca$cluster[ligneR[1]]
      ligneB <- tabClusterSca[which(tabClusterSca$cluster==clus),]
      ligneB <- ligneB$range
      ligneB <- ligneB[which(ligneB %in% ligneR)]
      ligneR <- ligneR[-which(ligneR %in% ligneB)]
      
      tabTemp <- tab[ligneB,]
      #find insertion point
      refPos <- c(NA, NA, NA, NA)
      if(!is.na(tabTemp$REF_POS[1])){ 
	refPos <- as.data.frame(table(tabTemp$REF_POS), stringsAsFactors=FALSE)
	lm <- localMaxima(refPos$Freq)
	if (length(lm)==1) {
	    refPos <- c(as.numeric(refPos[lm,1]), as.numeric(refPos[lm,2]), NA, NA)
	} else {
	    refPos <- refPos[lm,]    
	    refPos <- refPos[order(refPos$Freq,decreasing = TRUE),]
	    refPos <- c(as.numeric(refPos[1,1]), as.numeric(refPos[1,2]), as.numeric(refPos[2,1]), as.numeric(refPos[2,2]))
	}
	#refPos <- floor(sum(as.numeric(tab$REF_POS[ligneB]))/length(ligneB))
      }
      refStart <- min(c(tabTemp$REF_START, tabTemp$REF_END))
      refEnd <- max(c(tabTemp$REF_START, tabTemp$REF_END))

      scaStart <- min(c(tabTemp$SCA_START, tabTemp$SCA_END))
      scaEnd <- max(c(tabTemp$SCA_START, tabTemp$SCA_END))

      scaPos <- c(NA, NA, NA, NA)
      if(!is.na(tabTemp$SCA_POS[1])){
	scaPos <- as.data.frame(table(tabTemp$SCA_POS), stringsAsFactors=FALSE)
	lm <- localMaxima(scaPos$Freq)
	if (length(lm)==1) {
	    scaPos <- c(as.numeric(scaPos[lm,1]), as.numeric(scaPos[lm,2]), NA, NA)
	} else {
	    scaPos <- scaPos[lm,]    
	    scaPos <- scaPos[order(scaPos$Freq,decreasing = TRUE),]
	    scaPos <- c(as.numeric(scaPos[1,1]), as.numeric(scaPos[1,2]), as.numeric(scaPos[2,1]), as.numeric(scaPos[2,2]))
	}
	#scaPos <- floor(sum(as.numeric(tab$SCA_POS_A[ligneB]))/length(ligneB))
      } 

      twoStrand <- FALSE
      tabStrand <- tabTemp[!duplicated(tabTemp$REF_STRAND),]
      if(length(tabStrand[,1])>1){
	twoStrand <- TRUE
      }
      tabStrand <- tabTemp[!duplicated(tabTemp$SCA_STRAND),]
      if(length(tabStrand[,1])>1){
	twoStrand <- TRUE
      }

      cluster <- rbind(cluster, data.frame(CHR=tabTemp$CHR[1], REF_START=refStart, REF_END=refEnd, REF_POS_A=refPos[1], NB_REF_POS_A=refPos[2], REF_POS_B=refPos[3], NB_REF_POS_B=refPos[4], SCA_START=scaStart, SCA_END=scaEnd, SCA_POS_A=scaPos[1], NB_SCA_POS_A=scaPos[2], SCA_POS_B=scaPos[3], NB_SCA_POS_B=scaPos[4], SCA_READ_FIRST=first, NUM_READ=length(ligneB), TYPE=type, TWO_STRAND=twoStrand, stringsAsFactors=FALSE) )

    }

  }

  return (cluster)
}

selectOEAinsert <- function(file, param, idSca, type, first) {

  OEAsca <- readGAlignments(file, param=param,use.names=TRUE)
  seqlevels(OEAsca) <- as.character(idSca)

  idSC <- 1
  if(first) { idSC <- 2 }
  bamRef <- Gbam[[idSC]][which(names(Gbam[[idSC]]) %in% names(OEAsca))]

  bamRef <- bamRef[order(names(bamRef))]
  OEAsca <- OEAsca[order(names(OEAsca))]

# ##########################################################
# Récupère la partie du OEA qui est sur la reference
# ##########################################################

  if (length(OEAsca) > length(bamRef)) {
      read <- as.data.frame(table(names(OEAsca)), stringsAsFactors=FALSE)
      read2 <- read[which(read$Freq>1),]
      for (i in 1:length(read2[,1])){
	if (read2$Freq[i]==2) {
	  iRead <- grep(read2$Var1[i], names(OEAsca), fixed=TRUE)
	  badOne <- iRead[which(qwidth(OEAsca[iRead])==min(qwidth(OEAsca[iRead])))]
	  OEAsca <- OEAsca[-badOne]
	}else{
	  OEAsca <- OEAsca[-which(names(OEAsca)==read2$Var1[i])]
	  bamRef <- bamRef[-which(names(bamRef)==read2$Var1[i])]
	}
      }
  }
  if (length(OEAsca) < length(bamRef)) {
      read <- as.data.frame(table(names(bamRef)), stringsAsFactors=FALSE)
      read2 <- read[which(read$Freq>1),]
      for (i in 1:length(read2[,1])){
	if (read2$Freq[i]==2) {
	  iRead <- grep(read2$Var1[i], names(bamRef), fixed=TRUE)
	  badOne <- iRead[which(qwidth(bamRef[iRead])==min(qwidth(bamRef[iRead])))]
	  bamRef <- bamRef[-badOne]
	}else{
	  OEAsca <- OEAsca[-which(names(OEAsca)==read2$Var1[i])]
	  bamRef <- bamRef[-which(names(bamRef)==read2$Var1[i])]
	}
      }
  }

  if (length(OEAsca) > length(bamRef)) {
      read <- as.data.frame(table(names(OEAsca)), stringsAsFactors=FALSE)
      read2 <- read[which(read$Freq>1),]
      for (i in 1:length(read2[,1])){
	if (read2$Freq[i]==2) {
	  iRead <- grep(read2$Var1[i], names(OEAsca), fixed=TRUE)
	  badOne <- iRead[which(qwidth(OEAsca[iRead])==min(qwidth(OEAsca[iRead])))]
	  OEAsca <- OEAsca[-badOne]
	}else{
	  OEAsca <- OEAsca[-which(names(OEAsca)==read2$Var1[i])]
	  bamRef <- bamRef[-which(names(bamRef)==read2$Var1[i])]
	}
      }
  }

  if (length(OEAsca) != length(bamRef)) {# Shouldn't append but ...
    stop(paste("error 2 : number of OEAsca is different. Need to improve the R script. idSca :",idSca , sep=""))
  }

  if(length(OEAsca)==0){return(NULL)}

# #####################################
# clustering of overlapped read
# #####################################
  
  # merge map and unmap region of each OEA
  tabInsert <- data.frame(NAME=names(bamRef), CHR=as.vector(seqnames(bamRef)), REF_START=start(bamRef), REF_END=end(bamRef), REF_POS=rep(NA,length(end(bamRef))), REF_STRAND=as.vector(strand(bamRef)), SCA_ID=rep(idSca,length(end(OEAsca))), SCA_START=start(OEAsca), SCA_END=end(OEAsca), SCA_POS=rep(NA,length(end(OEAsca))), SCA_STRAND=as.vector(strand(OEAsca)), SCA_READ_FIRST=rep(first,length(end(OEAsca))), stringsAsFactors=FALSE) 

  rm(bamRef)
  rm(OEAsca)

  cluster <- creatInsertCluster(tabInsert, idSca, type, first)

  return(cluster)

}

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  x <- c(0,x,0)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y<-1
  return(y)
}

addInsert <- function (selectIns){

  tab <- data.frame(REF_POS=c(selectIns$REF_POS_A, selectIns$REF_POS_B),SCA_POS=c(selectIns$SCA_POS_A, selectIns$SCA_POS_B), NB=c(selectIns$NB_REF_POS_A, selectIns$NB_REF_POS_B), SCA_READ_FIRST=c(selectIns$SCA_READ_FIRST,selectIns$SCA_READ_FIRST), stringsAsFactors=FALSE)
  tab <- cbind(tab, data.frame(KEY=paste(tab$REF_POS,tab$SCA_POS,sep=":"), stringsAsFactors=FALSE))
  tab <- tab[-which(is.na(tab$REF_POS | tab$SCA_POS)),]

  allInsertTmp <- NULL

  if (!is.null(tab) & length(tab[,1])>0) {
    ins <- NULL
    for (l in unique(tab$KEY)) {
	selectL <- tab[which(tab$KEY==l),]
	ins <- rbind(ins,data.frame(REF_POS=selectL$REF_POS[1], SCA_POS=selectL$SCA_POS[1], NB=sum(selectL$NB), stringsAsFactors=FALSE))
    }
    rm(tab)

    ins <- ins[order(ins$NB, decreasing=TRUE),]
    if (length(ins[,1])<2) {
      refPos1 <- ins$REF_POS[1]
      scaPos1 <- ins$SCA_POS[1]
      nbInser1 <- ins$NB[1]
      refPos2 <- NA
      scaPos2 <- NA
      nbInser2 <- NA
    } else {
      refPos1 <- ins$REF_POS[1]
      scaPos1 <- ins$SCA_POS[1]
      nbInser1 <- ins$NB[1]
      refPos2 <- ins$REF_POS[2]
      scaPos2 <- ins$SCA_POS[2]
      nbInser2 <- ins$NB[2]
    }

    sclip1 <- sum(selectIns$NUM_READ[which(selectIns$TYPE=="sclip" & selectIns$SCA_READ_FIRST==TRUE)])
    sclip2 <- sum(selectIns$NUM_READ[which(selectIns$TYPE=="sclip" & selectIns$SCA_READ_FIRST==FALSE)])
    OEA1 <- sum(selectIns$NUM_READ[which(selectIns$TYPE=="OEA" & selectIns$SCA_READ_FIRST==TRUE)])
    OEA2 <- sum(selectIns$NUM_READ[which(selectIns$TYPE=="OEA" & selectIns$SCA_READ_FIRST==FALSE)])

    twoStrand <- 0
    if (sum(selectI$TWO_STRAND[which(selectI$TYPE=="sclip" & selectI$SCA_READ_FIRST==TRUE)])>=1) {
      twoStrand <- twoStrand + 1
    }
    if (sum(selectI$TWO_STRAND[which(selectI$TYPE=="sclip" & selectI$SCA_READ_FIRST==FALSE)])>=1) {
      twoStrand <- twoStrand + 1
    }
    if (sum(selectI$TWO_STRAND[which(selectI$TYPE=="OEA" & selectI$SCA_READ_FIRST==TRUE)])>=1) {
      twoStrand <- twoStrand + 1
    }
    if (sum(selectI$TWO_STRAND[which(selectI$TYPE=="OEA" & selectI$SCA_READ_FIRST==FALSE)])>=1) {
      twoStrand <- twoStrand + 1
    }

    allInsertTmp <- rbind(allInsertTmp, data.frame(CHR=selectI$CHR[1], REF_START=min(c(selectI$REF_START,selectI$REF_END)), REF_END=max(c(selectI$REF_START,selectI$REF_END)), INSERT_POS_A=refPos1, NB_INSERT_A=nbInser1, INSERT_POS_B=refPos2, NB_INSERT_B=nbInser2, SCLIP1=sclip1, SCOEA1=scOEA1, OEA1=OEA1, SCLIP2=sclip2, SCOEA2=scOEA2, OEA2=OEA2, ID_SCA=selectI$ID_SCA[1], SCA_START=min(c(selectI$SCA_START,selectI$SCA_END)), SCA_END=max(c(selectI$SCA_START,selectI$SCA_END)), INSERT_SCA_A=scaPos1, INSERT_SCA_B=scaPos2, TWO_STRAND=twoStrand, stringsAsFactors=FALSE ))
  }

  return(allInsertTmp)

}


clusterFusion <- function(tabInsert, insertSize) {
  l<-1
  while(l < length(tabInsert$CHR)){
    maxSca <- max(c(tabInsert$SCA_END[l+1], tabInsert$SCA_END[l], tabInsert$SCA_START[l], tabInsert$SCA_START[l+1]))
    minSca <- min(c(tabInsert$SCA_END[l+1], tabInsert$SCA_END[l], tabInsert$SCA_START[l], tabInsert$SCA_START[l+1]))
    maxRef <- max(c(tabInsert$REF_END[l+1], tabInsert$REF_END[l], tabInsert$REF_START[l], tabInsert$REF_START[l+1]))
    minRef <- min(c(tabInsert$REF_END[l+1], tabInsert$REF_END[l], tabInsert$REF_START[l], tabInsert$REF_START[l+1]))
    if ( (tabInsert$CHR[l]==tabInsert$CHR[l+1]) & (tabInsert$ID_SCA[l]==tabInsert$ID_SCA[l+1]) & (abs(maxRef-minRef)<(insertSize)) & (abs(maxSca-minSca)<(insertSize)) ) {
      tabInsert$REF_END[l] <- maxRef
      tabInsert$REF_START[l] <- minRef

      tabInsert$SCA_END[l] <- maxSca
      tabInsert$SCA_START[l] <- minSca

      tabInsert$NUM_READ[l] <- tabInsert$NUM_READ[l] + tabInsert$NUM_READ[l+1]

      tabInsert$TWO_STRAND[l] <- tabInsert$TWO_STRAND[l] | tabInsert$TWO_STRAND[l+1]

      tabInsert <- tabInsert[-(l+1),]
    } else {
      l<-l+1
    }
  }

  return(tabInsert)
}

# #########################################################################################################
# #########################################################################################################
# MAIN
# #########################################################################################################
# #########################################################################################################
folderOut

file <- paste(folderOut,"scaffolds.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}

tabSca <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character", "numeric", "numeric"))

allCluster <- NULL
idInsert <- 1

for (i in 1:length(tabSca$Scaffold)) {

  idSca <- tabSca$Scaffold[i]
  match <- tabSca$MATCH[i]
  param <- ScanBamParam(what="seq", which=GRanges(idSca, IRanges(1, tabSca$LEN_SCA[i])))

  seqInsertion <- c(1,tabSca$LEN_SCA[i])
  
  C_OEA1 <- selectOEAinsert(paste(folderScaffold,"cov/OEAUNMAP.1.bam",sep=""), param, idSca, "OEA", TRUE)
  C_OEA2 <- selectOEAinsert(paste(folderScaffold,"cov/OEAUNMAP.2.bam",sep=""), param, idSca, "OEA", FALSE)

  cluster <- rbind(C_OEA1, C_OEA2)
  
  if( length(cluster[,1])>0 ) {
    cluster <- cbind(cluster, data.frame(REGION=paste(idSca,"_",cluster$CHR, sep=""), ID_SCA=rep(idSca, length(cluster$CHR)), stringsAsFactors=FALSE))
    allCluster <- rbind(allCluster, cluster)
  }
}

if(!is.null(allCluster)){
  write.table(allCluster, file=paste(folderOut,"cluster.OEA.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)
  allCluster <- allCluster[order(allCluster$ID_SCA, allCluster$CHR, allCluster$REF_START), ]
  allClusterF <- allCluster[which(allCluster$SCA_READ_FIRST==FALSE),]
  allClusterT <- allCluster[which(allCluster$SCA_READ_FIRST==TRUE),]
  allClusterF <- clusterFusion(allClusterF, insertSize)
  allClusterT <- clusterFusion(allClusterT, insertSize)
  allCluster2 <- rbind(allClusterF, allClusterT)
  allCluster2 <- allCluster2[ order(allCluster2$ID_SCA, allCluster2$CHR, allCluster2$REF_START), ]
  write.table(allCluster2, file=paste(folderOut,"cluster.OEA.fusion.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)
}


