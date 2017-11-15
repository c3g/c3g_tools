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
insertSize <- as.numeric(args[5])
nbThread <- as.numeric(args[6])
maxPoint <- 500

algo <- ""
if (length(args)>6) {
  algo <- args[7]
}

#lib
library(Rsamtools)

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

flagReadMapped <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE, isUnmappedQuery=FALSE, hasUnmappedMate=FALSE, isDuplicate=FALSE)
flagReadMappedDisc <- scanBamFlag(isPaired=TRUE, isProperPair=FALSE, isUnmappedQuery=FALSE, hasUnmappedMate=FALSE, isDuplicate=FALSE)
flagReadOEAMap <- scanBamFlag(isPaired=TRUE, isUnmappedQuery=FALSE, hasUnmappedMate=TRUE, isDuplicate=FALSE)

# #########################################################################################################
# #########################################################################################################
# FONCTION
# #########################################################################################################
# #########################################################################################################


splitBadQuality <- function(bam){
  bamBadQuality <- NULL
  if(length(bam)>0){
    bamTmp <- as.data.frame(bam, row.names = as.character(c(1:length(bam))))
    bad <- which(bamTmp$mapq<10)
    if (length(bad)>0){
      bamBadQuality <- bam[bad,]
      bam <- bam[-bad,]
    }
  }
  return (list(bam, bamBadQuality))
}

computeCovRef <- function(bam, id, type, start, end) {
  cov <- coverage(bam)
  cov <- as.vector(cov[[as.character(seqnames(bam))[1]]])
  cov <- cov[start:end]
  covTab <- data.frame(ID=rep(id,length(cov)), POS=c(1:length(cov)), COV=cov, TYPE=rep(type,length(cov)), stringsAsFactors=FALSE)
  covTab$POS<-c(start:end)
  return(covTab)
}

covInsert <- function (j){
  start <- allInsert$REF_START[j]-insertSize*4
  end <- allInsert$REF_END[j]+insertSize*4

  if (!is.na(allInsert$INSERT_POS_A[j])) {
    pos <- c(allInsert$INSERT_POS_A[j], allInsert$INSERT_POS_B[j])
    pos <- pos[!is.na(pos)]
    start <- min(pos) - insertSize*6
    end <- max(pos) + insertSize*6
  }

  allCovIns <- NULL

  if (start<1) {start <- 1}
  paramReadMapped <- ScanBamParam(what=c("seq","mapq","isize"), flag=flagReadMapped, which=GRanges(allInsert$CHR[j], IRanges(start, end)))
  paramReadMappedDisc <- ScanBamParam(what=c("seq","mapq","isize"), flag=flagReadMappedDisc, which=GRanges(allInsert$CHR[j], IRanges(start, end)))
  paramReadOEAMap <- ScanBamParam(what=c("seq","mapq","isize"), flag=flagReadOEAMap, which=GRanges(allInsert$CHR[j], IRanges(start, end)))

  bamMapped <- readGAlignments(paste(folderSC,sample,".scOthers.bam",sep=""),param=paramReadMapped,use.names=TRUE)
  bamOEAMap <- readGAlignments(paste(folderExtract,sample,".OEAMAP.bam",sep=""),param=paramReadOEAMap,use.names=TRUE)
  bamScOEAMap <- readGAlignments(paste(folderExtract,sample,".scOEAMAP.bam",sep=""),param=paramReadOEAMap,use.names=TRUE)
  bamDiscordant <- readGAlignments(paste(folderSC,sample,".scOthers.bam",sep=""),param=paramReadMappedDisc,use.names=TRUE)
  bamMappedBadQuality <- NULL
  bamOEAMapBadQuality <- NULL
  bamScOEAMapBadQuality <- NULL

  resList <- splitBadQuality(bamMapped)
  bamMapped <- resList[[1]]
  bamMappedBadQuality <- resList[[2]]

  resList <- splitBadQuality(bamOEAMap)
  bamOEAMap <- resList[[1]]
  bamOEAMapBadQuality <- resList[[2]]

  resList <- splitBadQuality(bamScOEAMap)
  bamScOEAMap <- resList[[1]]
  bamScOEAMapBadQuality <- resList[[2]]

  idInsert <- allInsert$ID_INSERT[j]

  if (length(bamMapped)>0){
    cov <- computeCovRef(bamMapped, idInsert, "Mapped", start, end)
    allCovIns <- rbind(allCovIns, cov)
  }
  if (length(bamDiscordant)>0){
    cov <- computeCovRef(bamDiscordant, idInsert, "Discordant", start, end)
    allCovIns <- rbind(allCovIns, cov)
  }
  if (length(bamOEAMap)>0){
    cov <- computeCovRef(bamOEAMap, idInsert, "OEA", start, end)
    allCovIns <- rbind(allCovIns, cov)
  }
  if (length(bamScOEAMap)>0){
    cov <- computeCovRef(bamScOEAMap, idInsert, "scOEA", start, end)
    allCovIns <- rbind(allCovIns, cov)
  }

  if (!is.null(bamMappedBadQuality)){
    cov <- computeCovRef(bamMappedBadQuality, idInsert, "BadQ", start, end)
    allCovIns <- rbind(allCovIns, cov)
  }
  if (!is.null(bamOEAMapBadQuality)){
    cov <- computeCovRef(bamOEAMapBadQuality, idInsert, "OEA BQ", start, end)
    allCovIns <- rbind(allCovIns, cov)
  }
  if (!is.null(bamScOEAMapBadQuality)){
    cov <- computeCovRef(bamScOEAMapBadQuality, idInsert, "scOEA BQ", start, end)
    allCovIns <- rbind(allCovIns, cov)
  }
  return(allCovIns)
}

# #########################################################################################################
# #########################################################################################################
# MAIN
# #########################################################################################################
# #########################################################################################################

file <- paste(folderOut,"insert.filter.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}

allInsert <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

if (length(allInsert$CHR)==0){
q()
}

param <- ScanBamParam(what="seq")

allCovIns <- data.frame(ID_INS="1", POS=1, COV=1, TYPE="Mapped")
allCovIns <- allCovIns[-1,]
write.table(allCovIns, file=paste(folderOut,"covIns.filter.tab",sep=""), sep="\t", dec = ".", row.names=FALSE, col.names=TRUE, quote=FALSE)

allCovIns <- data.frame(ID_INS="1", POS=1, COV=1, MIN=1, MAX=1, TYPE="Mapped")
allCovIns <- allCovIns[-1,]
write.table(allCovIns, file=paste(folderOut,"covIns.filter",maxPoint,".tab",sep=""), sep="\t", dec = ".", row.names=FALSE, col.names=TRUE, quote=FALSE)

i <- 1
while( i <= length(allInsert$CHR)){
  #create table of covIns
  max <- i+nbThread-1
  if(max>length(allInsert$CHR))	{max <- length(allInsert$CHR)}
  tabIndex <- c(i:max)
  for(j in tabIndex){
    cmd=expression(covInsert(j))
    mcparallel(cmd)
  }
  tmp <- mccollect()

  for (j in c(1:length(tabIndex))){
    allCovIns <- tmp[[j]]
    if (!is.null(allCovIns)) {
      colnames(allCovIns)[1] <- "ID_INS"
      write.table(allCovIns, file=paste(folderOut,"covIns.filter.tab",sep=""), sep="\t", dec=".", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
      meanCov <- tapply(allCovIns$COV,cut(1:length(allCovIns$TYPE),maxPoint),mean)
      minCov <- tapply(allCovIns$COV,cut(1:length(allCovIns$TYPE),maxPoint),min)
      maxCov <- tapply(allCovIns$COV,cut(1:length(allCovIns$TYPE),maxPoint),max)
      allCovIns <- data.frame(ID=rep(allCovIns$ID[1],maxPoint), POS=round(seq(allCovIns$POS[1], allCovIns$POS[length(allCovIns$POS)],length.out=maxPoint),0), COV=meanCov, MIN=minCov, MAX=maxCov, TYPE=rep(allCovIns$TYPE[1],maxPoint), stringsAsFactors=FALSE)
      write.table(allCovIns, file=paste(folderOut,"covIns.filter",maxPoint,".tab",sep=""), sep="\t", dec=".", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    }
  }

  i <- max+1
}

#write.table(allCovIns, file=paste(folderOut,"covIns.filter.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)





