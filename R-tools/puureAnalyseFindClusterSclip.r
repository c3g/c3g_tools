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
minMapq <- as.numeric(args[5])
insertSize <- as.numeric(args[6])

algo <- ""
if (length(args)>6) {
  algo <- args[7]
}

#lib
library(Rsamtools)
library(grid)
library(igraph) # install.packages("igraph", lib="~/.R/library")

#set path 
folderScaffold <- paste(mainPathSca,"/scaffolds",algo,"/",sample,"/ray/ray21/",sep="")
folderSC <- paste(mainPath,"/sclip",algo,"/",sample,"/",sep="")
folderExtract <- paste(mainPath,"/extract",algo,"/",sample,"/",sep="")
folderOut <- paste(folderScaffold,"/insert",type,"/", sep="")
dir.create(file.path(folderOut), showWarnings = FALSE, recursive=TRUE)


# #########################################################################################################
# #########################################################################################################
# GLOBAL
# #########################################################################################################
# #########################################################################################################
flagFirstRead<- scanBamFlag(isFirstMateRead=TRUE)
flagSecondRead<- scanBamFlag(isFirstMateRead=FALSE)

paramFirst <- ScanBamParam(what="seq", flag=flagFirstRead)
paramSecond <- ScanBamParam(what="seq", flag=flagSecondRead)

Gbam <- list(readGAlignments(paste(folderSC,sample,".sc.bam",sep=""),param=paramFirst,use.names=TRUE), readGAlignments(paste(folderSC,sample,".sc.bam",sep=""),param=paramSecond,use.names=TRUE))

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

selectSclipInsert <- function(file, param, idSca, type, first) {

  sclipSca <- readGAlignments(file, param=param,use.names=TRUE)
  seqlevels(sclipSca) <- as.character(idSca)

  idSC <- 1
  if(!first) { idSC <- 2 }
  bamSclip <- Gbam[[idSC]][which(names(Gbam[[idSC]]) %in% names(sclipSca))]
  bamSclip <- bamSclip[grep("S", cigar(bamSclip), fixed=TRUE)]

  bamSclip <- bamSclip[order(names(bamSclip))]
  sclipSca <- sclipSca[order(names(sclipSca))]

# ##########################################################
# Récupère la partie du softclip qui est sur la reference
# ##########################################################

  if (length(sclipSca) > length(bamSclip)) {
      read <- as.data.frame(table(names(sclipSca)), stringsAsFactors=FALSE)
      read2 <- read[which(read$Freq>1),]
      for (i in 1:length(read2[,1])) {
	iRead <- grep(read2$Var1[i], names(bamSclip), fixed=TRUE)
	jRead <- grep(read2$Var1[i], names(sclipSca), fixed=TRUE)
	bamSclip <- bamSclip[-iRead]
	sclipSca <- sclipSca[-jRead]
      }
  }

  if (length(sclipSca) < length(bamSclip)) {
      read <- as.data.frame(table(names(bamSclip)), stringsAsFactors=FALSE)
      read2 <- read[which(read$Freq==2),]
      for (i in 1:length(read2[,1])) {
	iRead <- grep(read2$Var1[i], names(bamSclip), fixed=TRUE)
	jRead <- grep(read2$Var1[i], names(sclipSca), fixed=TRUE)
	if(length(jRead)==1){
	  sclipScaTmp <- as.data.frame(sclipSca, row.names = as.character(c(1:length(sclipSca))))
	  bamSclipTmp <- as.data.frame(bamSclip, row.names = as.character(c(1:length(bamSclip))))
	  badOne <- iRead[grep(sclipScaTmp$seq[jRead], bamSclipTmp$seq[iRead], fixed=TRUE, invert=TRUE)]
	  if(length(badOne)==0){
	    stop(paste("error : 2 reads paired have same subseq softcliped. Need to improve the R script. idSca :",idSca , sep=""))
	  }
	  bamSclip <- bamSclip[-badOne]
	} else {
	  bamSclip <- bamSclip[-iRead]
	  sclipSca <- sclipSca[-jRead]
	}
      }
  }

  if (length(sclipSca) > length(bamSclip)) {
      read <- as.data.frame(table(names(sclipSca)), stringsAsFactors=FALSE)
      read2 <- read[which(read$Freq>1),]
      for (i in 1:length(read2[,1])) {
	iRead <- grep(read2$Var1[i], names(bamSclip), fixed=TRUE)
	jRead <- grep(read2$Var1[i], names(sclipSca), fixed=TRUE)
	bamSclip <- bamSclip[-iRead]
	sclipSca <- sclipSca[-jRead]
      }
  }


  if (length(sclipSca) != length(bamSclip)) { # it should be not append but ...
    stop(paste("error : 2 reads paired have softclip. Need to improve the R script. idSca :",idSca , sep=""))
  }

# #####################################
# clustering of overlapped read
# #####################################
  
  # merge map and unmap region of each softclip
  tabSclip <- data.frame(NAME=names(bamSclip), CHR=as.vector(seqnames(bamSclip)), REF_START=start(bamSclip), REF_END=end(bamSclip), REF_POS=end(bamSclip), REF_STRAND=as.vector(strand(bamSclip)), CIGAR=substr(cigar(bamSclip), nchar(cigar(bamSclip)), nchar(cigar(bamSclip))), SCA_ID=rep(idSca,length(end(sclipSca))), SCA_START=start(sclipSca), SCA_END=end(sclipSca), SCA_POS=start(sclipSca), SCA_STRAND=as.vector(strand(sclipSca)), SCA_READ_FIRST_READ=rep(first,length(end(sclipSca))), INSERT=paste(as.vector(seqnames(bamSclip)),end(bamSclip),sep=":"), stringsAsFactors=FALSE) 

  # select position os sclip (end or start of the alignment)
  i <-which(tabSclip$CIGAR=="M")
  tabSclip$REF_POS[i] <- tabSclip$REF_START[i]
  tabSclip$SCA_POS[i] <- tabSclip$SCA_END[i]
  tabSclip$INSERT[i] <- paste(as.vector(seqnames(bamSclip))[i],tabSclip$REF_START[i],sep=":")

  rm(bamSclip)
  rm(sclipSca)

  cluster <- creatInsertCluster(tabSclip, idSca, type, first)

  tabSclip <- tabSclip[c(1,2,5,8,12,13),]

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
  
  C_SC1 <- selectSclipInsert(paste(folderScaffold,"cov/sclip.1.bam",sep=""), param, idSca, "sclip", TRUE)
  C_SC2 <- selectSclipInsert(paste(folderScaffold,"cov/sclip.2.bam",sep=""), param, idSca, "sclip", FALSE)
  cluster <- rbind(C_SC1, C_SC2)

# no efficient
#   cmd=expression(selectSclipInsert(paste(folderScaffold,"cov/sclip.1.bam",sep=""), param, idSca, "sclip", TRUE))
#   mcparallel(cmd)
#   cmd=expression(selectSclipInsert(paste(folderScaffold,"cov/sclip.2.bam",sep=""), param, idSca, "sclip", FALSE))
#   mcparallel(cmd)
#   tmp <- mccollect()
#   cluster <- rbind(tmp[[1]], tmp[[2]])
  
  if( length(cluster[,1])>0 ) {
    cluster <- cbind(cluster, data.frame(REGION=paste(idSca,"_",cluster$CHR, sep=""), ID_SCA=rep(idSca, length(cluster$CHR)), stringsAsFactors=FALSE))
    allCluster <- rbind(allCluster, cluster)
  }
}

if (!is.null(allCluster)) {
  write.table(allCluster, file=paste(folderOut,"cluster.sc.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)
  allCluster <- allCluster[order(allCluster$ID_SCA, allCluster$CHR, allCluster$REF_START), ]
  allClusterF <- allCluster[which(allCluster$SCA_READ_FIRST==FALSE),]
  allClusterT <- allCluster[which(allCluster$SCA_READ_FIRST==TRUE),]
  allClusterF <- clusterFusion(allClusterF, insertSize)
  allClusterT <- clusterFusion(allClusterT, insertSize)
  allCluster2 <- rbind(allClusterF, allClusterT)
  allCluster2 <- allCluster2[ order(allCluster2$ID_SCA, allCluster2$CHR, allCluster2$REF_START), ]
  write.table(allCluster2, file=paste(folderOut,"cluster.sc.fusion.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)
}




