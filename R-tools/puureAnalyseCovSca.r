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
nbThread <- as.numeric(args[5])
maxPoint <- 500
algo <- ""
if (length(args)>5) {
  algo <- args[6]
}

#lib
library(Rsamtools)

#set path 
folderScaffold <- paste(mainPathSca,"/scaffolds",algo,"/",sample,"/ray/ray21/",sep="")
folderOut <- paste(folderScaffold,"/insert",type,"/", sep="")
dir.create(file.path(folderOut), showWarnings = FALSE, recursive=TRUE)

# #########################################################################################################
# #########################################################################################################
# FONCTION
# #########################################################################################################
# #########################################################################################################

computeCovSca <- function(bam, id) {
  cov <- coverage(bam)
  cov <- as.vector(cov[[as.character(id)]])
  return(cov)
}

computeAllCov <- function(folderScaffold, param, idSca) {

  sclip1 <- readGAlignments(paste(folderScaffold,"cov/sclip.1.bam",sep=""),param=param,use.names=TRUE)
  sclip2 <- readGAlignments(paste(folderScaffold,"cov/sclip.2.bam",sep=""),param=param,use.names=TRUE)
  scOEA1 <- readGAlignments(paste(folderScaffold,"cov/scOEAUNMAP.1.bam",sep=""),param=param,use.names=TRUE)
  scOEA2 <- readGAlignments(paste(folderScaffold,"cov/scOEAUNMAP.2.bam",sep=""),param=param,use.names=TRUE)
  OEA1 <- readGAlignments(paste(folderScaffold,"cov/OEAUNMAP.1.bam",sep=""),param=param,use.names=TRUE)
  OEA2 <- readGAlignments(paste(folderScaffold,"cov/OEAUNMAP.2.bam",sep=""),param=param,use.names=TRUE)
  orphan <- readGAlignments(paste(folderScaffold,"cov/ORPHAN.bam",sep=""),param=param,use.names=TRUE)
  seqlevels(OEA1) = idSca
  seqlevels(OEA2) = idSca
  seqlevels(scOEA1) = idSca
  seqlevels(scOEA2) = idSca
  seqlevels(sclip1) = idSca
  seqlevels(sclip2) = idSca
  seqlevels(orphan) = idSca

  cov <- computeCovSca(orphan, idSca)
  covTab <- data.frame(ID=rep(idSca,length(cov)), POS=c(1:length(cov)), COV=cov, TYPE=rep("Orphan",length(cov)), stringsAsFactors=FALSE)
  
  cov1 <- computeCovSca(sclip1, idSca)  
  cov2 <- computeCovSca(sclip2, idSca)
  cov <- cov1 + cov2
  covTab <- rbind(covTab,data.frame(ID=rep(idSca,length(cov)), POS=c(1:length(cov)), COV=cov, TYPE=rep("Sclip",length(cov)), stringsAsFactors=FALSE))

  cov1 <- computeCovSca(scOEA1, idSca)
  cov2 <- computeCovSca(scOEA2, idSca)
  cov <- cov1 + cov2
  covTab <- rbind(covTab,data.frame(ID=rep(idSca,length(cov)), POS=c(1:length(cov)), COV=cov, TYPE=rep("scOEA",length(cov)), stringsAsFactors=FALSE))
  
  cov1 <- computeCovSca(OEA1, idSca)
  cov2 <- computeCovSca(OEA2, idSca)
  cov <- cov1 + cov2
  covTab <- rbind(covTab,data.frame(ID=rep(idSca,length(cov)), POS=c(1:length(cov)), COV=cov, TYPE=rep("OEA",length(cov)), stringsAsFactors=FALSE))
  
  return(covTab)
}

covSca <- function (i) {

  idSca <- tabSca$ID_SCA[i]
  match <- tabSca$MATCH[i]
  param <- ScanBamParam(what="seq", which=GRanges(idSca, IRanges(1, tabSca$LEN_SCA[i])))

  allCov <- computeAllCov(folderScaffold, param, idSca)
  return(allCov)
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

tabSca <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character", "numeric", "numeric"))

if (length(tabSca$ID_SCA)==0){
q()
}


param <- ScanBamParam(what="seq")


allCovSca <- data.frame(ID_SCA="1", POS=1, COV=1, TYPE="Init")
allCovSca <- allCovSca[-1,]
write.table(allCovSca, file=paste(folderOut,"covSca.filter.tab",sep=""), sep="\t", dec = ".", row.names=FALSE, col.names=TRUE, quote=FALSE)
allCovSca <- data.frame(ID_SCA="1", POS=1, COV=1, MIN=1, Max=1, TYPE="Init")
allCovSca <- allCovSca[-1,]
write.table(allCovSca, file=paste(folderOut,"covSca.filter",maxPoint,".tab",sep=""), sep="\t", dec = ".", row.names=FALSE, col.names=TRUE, quote=FALSE)

i <- 1
while( i <= length(tabSca$ID_SCA)){
  #create table of covIns
  max <- i+nbThread-1
  if(max>length(tabSca$ID_SCA))	{max <- length(tabSca$ID_SCA)}
  tabIndex <- c(i:max)
  for(j in tabIndex){
    cmd=expression(covSca(j))
    mcparallel(cmd)
  }
  tmp <- mccollect()

  for (j in c(1:length(tabIndex))){
    allCovSca <- tmp[[j]]
    if(!is.null(allCovSca)){
      colnames(allCovSca)[1] <- "ID_SCA"
      write.table(allCovSca, file=paste(folderOut,"covSca.filter.tab",sep=""), sep="\t", dec=".", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
      meanCov <- tapply(allCovSca$COV,cut(1:length(allCovSca$TYPE),maxPoint),mean)
      minCov <- tapply(allCovSca$COV,cut(1:length(allCovSca$TYPE),maxPoint),min)
      maxCov <- tapply(allCovSca$COV,cut(1:length(allCovSca$TYPE),maxPoint),max)
      allCovSca <- data.frame(ID=rep(allCovSca$ID[1],maxPoint), POS=round(seq(allCovSca$POS[1], allCovSca$POS[length(allCovSca$POS)],length.out=maxPoint),0), COV=meanCov, MIN=minCov, MAX=maxCov, TYPE=rep(allCovSca$TYPE[1],maxPoint), stringsAsFactors=FALSE)
      write.table(allCovSca, file=paste(folderOut,"covSca.filter",maxPoint,".tab",sep=""), sep="\t", dec=".", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    }
  }

  i <- max+1
}


