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
meanCov <- as.numeric(args[5])
insertSize <- as.numeric(args[6])
minOverlap <- as.numeric(args[7])
fileExclu <- args[8]

algo <- ""
if (length(args)>8) {
  algo <- args[9]
}

#lib
library(Rsamtools)
library(igraph) # install.packages("igraph", lib="~/.R/library")

#set path 
folderScaffold <- paste(mainPathSca,"/scaffolds",algo,"/",sample,"/ray/ray21/",sep="")
folderOut <- paste(folderScaffold,"/insert",type,"/", sep="")
dir.create(file.path(folderOut), showWarnings = FALSE, recursive=TRUE)

# #########################################################################################################
# #########################################################################################################
# FONCTION
# #########################################################################################################
# #########################################################################################################

addScoreMap <- function(tab, exclRegions, meanCov, insertSize){

  score <- tab$SCLIP1 + tab$SCOEA1 + tab$SCLIP2 + tab$SCOEA2
  tab <- cbind(tab, data.frame(SCORE=score, stringsAsFactors=FALSE))

  tmpInExclu <- tab[which(tab$CHR %in% unique(exclRegions[,1])),]
  tmpNotInExclu <- tab[which(!(tab$CHR %in% unique(exclRegions[,1]))),]
  if(length(tmpInExclu[,1])>0) {
    tmpInExclu <- cbind(tmpInExclu, data.frame(MAP=rep(1,length(tmpInExclu[,1]))))
  }
  if(length(tmpNotInExclu[,1])>0) {
    tmpNotInExclu <- cbind(tmpNotInExclu, data.frame(MAP=rep(NA,length(tmpNotInExclu[,1]))))
  }

  gr = GRanges(tmpInExclu$CHR,IRanges(start=(tmpInExclu$REF_START-insertSize),end=(tmpInExclu$REF_END+insertSize)))
  er.gr = GRanges(exclRegions[,1],IRanges(exclRegions[,2],exclRegions[,3]))

  fo = findOverlaps(gr,er.gr)
  if(length(fo)!=0){
    covOl = width(ranges(pintersect(gr[queryHits(fo)],er.gr[subjectHits(fo)]))) / width(gr[queryHits(fo)])
    ol.df = data.frame(index=queryHits(fo),cov.Ol=covOl)
    cov.ex = aggregate(covOl ~ index, data=ol.df,sum)
    if(length(tmpInExclu[,1])>0) {
      tmpInExclu$mappability[cov.ex$index] <- tmpInExclu$mappability[cov.ex$index] - cov.ex$covOl
    }
  }

  tab <- rbind(tmpInExclu,tmpNotInExclu)


  return(tab)
}

addInsert <- function (selectIns, minOverlap){

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
    ins <- ins[order(ins$NB, decreasing=TRUE),]
    if (length(ins[,1])<minOverlap) {
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
    scOEA1 <- sum(selectIns$NUM_READ[which(selectIns$TYPE=="scOEA" & selectIns$SCA_READ_FIRST==TRUE)])
    scOEA2 <- sum(selectIns$NUM_READ[which(selectIns$TYPE=="scOEA" & selectIns$SCA_READ_FIRST==FALSE)])

    twoStrand <- 0
    if (sum(selectI$TWO_STRAND[which(selectI$TYPE=="sclip" & selectI$SCA_READ_FIRST==TRUE)])>=1) {
      twoStrand <- twoStrand + 1
    }
    if (sum(selectI$TWO_STRAND[which(selectI$TYPE=="sclip" & selectI$SCA_READ_FIRST==FALSE)])>=1) {
      twoStrand <- twoStrand + 1
    }
    if (sum(selectI$TWO_STRAND[which(selectI$TYPE=="scOEA" & selectI$SCA_READ_FIRST==FALSE)])>=1) {
      twoStrand <- twoStrand + 1
    }
    if (sum(selectI$TWO_STRAND[which(selectI$TYPE=="scOEA" & selectI$SCA_READ_FIRST==TRUE)])>=1) {
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
    if ( (tabInsert$CHR[l]==tabInsert$CHR[l+1]) & (abs(tabInsert$REF_END[l+1]-tabInsert$REF_START[l])<(insertSize*2)) ) {
      tabInsert$REF_END[l] <- tabInsert$REF_END[l+1]
      tabInsert$INSERT_POS_B[l] <- tabInsert$INSERT_POS_A[l+1]
      tabInsert$NB_INSERT_B[l] <- tabInsert$NB_INSERT_A[l+1]

      tabInsert$SCLIP1[l] <- tabInsert$SCLIP1[l] + tabInsert$SCLIP1[l+1]
      tabInsert$SCOEA1[l] <- tabInsert$SCOEA1[l] + tabInsert$SCOEA1[l+1]
      tabInsert$OEA1[l] <- tabInsert$OEA1[l] + tabInsert$OEA1[l+1]
      tabInsert$SCLIP2[l] <- tabInsert$SCLIP2[l] + tabInsert$SCLIP2[l+1]
      tabInsert$SCOEA2[l] <- tabInsert$SCOEA2[l] + tabInsert$SCOEA2[l+1]
      tabInsert$OEA2[l] <- tabInsert$OEA2[l] + tabInsert$OEA2[l+1]

      tabInsert$SCA_END[l] <- tabInsert$SCA_END[l+1]
      tabInsert$INSERT_SCA_B[l] <- tabInsert$INSERT_SCA_A[l+1]
      tabInsert$TWO_STRAND[l] <- tabInsert$TWO_STRAND[l] + tabInsert$TWO_STRAND[l+1]
      tabInsert <- tabInsert[-(l+1),]
    }

    l<-l+1
  }

  return(tabInsert)
}


# #########################################################################################################
# #########################################################################################################
# MAIN
# #########################################################################################################
# #########################################################################################################

tabCluster<-NULL

file <- paste(folderOut,"cluster.sc.fusion.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
tabSC <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "logical", "numeric", "character", "logical", "character", "character"))

file <- paste(folderOut,"cluster.scOEA.fusion.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
tabSCOEA <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "logical", "numeric", "character", "logical", "character", "character"))

file <- paste(folderOut,"cluster.OEA.fusion.tab",sep="")
if (file.exists(file) & file.info(file)$size!=0){
  tabOEA <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "logical", "numeric", "character", "logical", "character", "character"))
  tabCluster <- rbind(tabSC, rbind(tabSCOEA, tabOEA))
} else {
  tabCluster <- rbind(tabSC, tabSCOEA)
}


rm(tabSC)
rm(tabSCOEA)
rm(tabOEA)

allInsert <- NULL
allCluster <- NULL
idInsert <- 1

allRegion <- unique(tabCluster$REGION)

for (region in allRegion) {
  cluster <- tabCluster[which(tabCluster$REGION==region),]
  if( length(cluster[,1])>0 ) {
    gr <-GRanges(seqnames=cluster$CHR, ranges = IRanges(start=cluster$REF_START, end=cluster$REF_END))
    g<-graph.edgelist(as.matrix(findOverlaps(gr)), directed=TRUE)
    insertRegion <- data.frame(cluster=as.numeric(clusters(g)$membership), range=as.numeric(V(g)), key=paste(as.numeric(clusters(g)$membership),as.numeric(V(g)),sep="-"), stringsAsFactors=FALSE)
    cpt <- 1
    for( j in unique(insertRegion[,1]) ){
      ligne <- insertRegion[which(insertRegion[,1]==j),]
      ligne <- ligne[,2]
      chr <- cluster$CHR[ligne[1]]
      refStart <- min(c(cluster$REF_START[ligne], cluster$REF_END[ligne]))
      refEnd <- max(c(cluster$REF_START[ligne], cluster$REF_END[ligne]))
      cluster$REGION[ligne] <- paste(chr,":", refStart, "-", refEnd, sep="")
      cpt <- cpt + 1
    }

    allInsertTmp <- NULL
    for ( j in unique(cluster$REGION) ){
      selectI <- cluster[which(cluster$REGION==j),]
      allInsertTmp <- rbind(allInsertTmp, addInsert(selectI, minOverlap))
    }

    if(!is.null(allInsertTmp)){
      allInsertTmp <- clusterFusion(allInsertTmp, insertSize)
      allInsert <- rbind(  allInsert, cbind( data.frame(ID_INSERT=seq(idInsert, (idInsert+length(allInsertTmp$CHR)-1)), stringsAsFactors=FALSE), allInsertTmp )  )

      cluster <- cbind( data.frame(ID_INSERT=rep(-1,length(cluster$CHR)), stringsAsFactors=FALSE),cluster )
      for ( j in 1:length(allInsertTmp$CHR)){
	cluster$ID_INSERT[which(cluster$CHR==allInsertTmp$CHR[j] & cluster$REF_START>=allInsertTmp$REF_START[j] & cluster$REF_START<=allInsertTmp$REF_END[j])] <- idInsert
	idInsert <- idInsert + 1
      }

      allCluster <- rbind(allCluster, cluster)
    }
  } 
}


if (is.null(allInsert)){
q()
}

allInsert <- allInsert[,c(1,15,2,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20)]
allCluster <- allCluster[,c(1,20,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]


exclRegions <- read.table(fileExclu)
if (!is.numeric(exclRegions[1,2])) {
exclRegions <- read.table(fileExclu,header=TRUE)
}

allInsert <- addScoreMap(allInsert, exclRegions, meanCov, insertSize)

allCluster <- allCluster[order(allCluster$ID_INSERT),]

write.table(allInsert, file=paste(folderOut,"insert.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)
write.table(allCluster, file=paste(folderOut,"cluster.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)


