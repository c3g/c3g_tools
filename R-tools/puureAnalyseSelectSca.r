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
minSize <- as.numeric(args[5])

algo <- ""
if (length(args)>5) {
  algo <- args[6]
}


#lib
library(IRanges)

#set path 
folderScaffold <- paste(mainPathSca,"/scaffolds",algo,"/",sample,"/ray/ray21/",sep="")
folderOut <- paste(folderScaffold,"/insert",type,"/", sep="")
dir.create(file.path(folderOut), showWarnings = FALSE, recursive=TRUE)


# #########################################################################################################
# #########################################################################################################
# MAIN
# #########################################################################################################
# #########################################################################################################

file <- paste(folderScaffold,"cov/readunmap.cov.txt", sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
cov <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

file <- paste(folderScaffold,"Scaffolds.fasta.refGenome.noHeader.psl" ,sep="")
if (!file.exists(file)) {
q()
}
blatDel <- NULL
if (file.info(file)$size!=0){
  tabBlat <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=FALSE, colClasses=c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character", "numeric", "numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "character", "character", "character")) 
  idSca <- unique(tabBlat[,10])

  blat <-NULL
  for (id in idSca) {
    blatTmp <- tabBlat[which(tabBlat[,10]==id),]
    blatTmp <- blatTmp[which(abs(blatTmp[,13]-blatTmp[,12])==max(abs(blatTmp[,13]-blatTmp[,12]))),]
    rangeTmp <- IRanges(blatTmp[,12], blatTmp[,13])
    rangeTmp <- as.data.frame(reduce(rangeTmp))
    nbCov <- sum(rangeTmp$width) - sum(blatTmp[,6])
    overRef <- data.frame(ID=id, REF=(nbCov/blatTmp[1,11])*100, stringsAsFactors=FALSE)
    blat <- rbind(blat, overRef)
  }
  blatDel <- blat[which(blat$REF>90),]
}


file <- paste(folderScaffold,"Scaffolds.fasta.blast",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
blast <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "numeric", "character", "character", "numeric", "numeric"))
blast <- blast[which(blast$T_MATCH==type),]


if(length(blast$ID)==0){
q()
}

# o <- which(is.na(cov$MeanCoverage))
# if (length(o)!=0){
#   cov <- cov[-o,]
# }
# 
# o <- which(cov$MeanCoverage<minCov)
# if (length(o)!=0){
#   cov <- cov[-o,]
# }
# 
# o <- which(blast$EVALUE>minEvalue)
# if (length(o)!=0){
#   blast <- blast[-o,]
# }


tabSca <-  merge(cov, blast, by.x="IntervalName", by.y="ID")
tabSca <- tabSca [,-c(2,3)]
colnames(tabSca) <- c("Scaffold", colnames(tabSca)[2:length(tabSca[1,])])

if(length(tabSca$Scaffold)==0){
q()
}

tabSca <- tabSca[order(tabSca$L_ID,decreasing=TRUE),]
colnames(tabSca)[14] <- "LEN_SCA"


o <- which(is.na(tabSca$MeanCoverage))
if (length(o)!=0){
  tabSca <- tabSca[-o,]
}
write.table(tabSca, file=paste(folderOut,"scaffolds.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)


nbRow <- length(tabSca$Scaffold)
#tabFil <- data.frame(ID_SCA=tabSca$Scaffold, COVREF=rep(TRUE,nbRow), MINCOV=rep(TRUE,nbRow), PROPERPAIRS=rep(TRUE,nbRow), MINLEN=rep(TRUE,nbRow), stringsAsFactors=FALSE)
tabFil <- data.frame(ID_SCA=tabSca$Scaffold, COVREF=rep(TRUE,nbRow), PROPERPAIRS=rep(TRUE,nbRow), MINLEN=rep(TRUE,nbRow), stringsAsFactors=FALSE)

#Del scaffolds who match on Ref
if (!is.null(blatDel)) {
  o <- which(tabSca$Scaffold %in% blatDel$ID)
  tabFil$COVREF[o] <- FALSE
}

#Del scaffolds with meanCov < minCov
#o <- which(tabSca$MeanCoverage<minCov)
#tabFil$MINCOV[o] <- FALSE

#Del scaffolds with not enough ProperPair
o <- which(tabSca$LEN_SCA>minSize & tabSca$NbPairedReads>0 & tabSca$NbProperPairs/tabSca$NbPairedReads<=0.50)
if (length(o)>0) {
  tabFil$PROPERPAIRS[o] <- FALSE
}
o <- which(tabSca$LEN_SCA>minSize & tabSca$NbPairedReads==0)
if (length(o)>0) {
  tabFil$PROPERPAIRS[o] <- FALSE
}

#Del scaffolds to small
o <- which(tabSca$LEN_SCA<=minSize)
tabFil$MINLEN[o] <- FALSE

total <- length(tabFil$ID_SCA)
nbCovRef <- length(which(tabFil$COVREF==FALSE))
nbMinCov <- length(which(tabFil$MINCOV==FALSE))
nbProperPairs <- length(which(tabFil$PROPERPAIRS==FALSE))
nbMinLen <- length(which(tabFil$MINLEN==FALSE))
good <- length(which(tabFil$COVREF==TRUE & tabFil$MINCOV==TRUE & tabFil$PROPERPAIRS==TRUE))
veryGood <- length(which(tabFil$COVREF==TRUE & tabFil$MINCOV==TRUE & tabFil$PROPERPAIRS==TRUE & tabFil$MINLEN==TRUE ))

total
nbCovRef
nbMinCov
nbProperPairs
nbMinLen
good
veryGood

write.table(tabFil, file=paste(folderOut,"scaffolds.toDelete.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)


# selectSca <- tabSca[which(tabSca$Scaffold %in% tabFil$ID_SCA[which(tabFil$COVNA==TRUE & tabFil$COVREF==TRUE & tabFil$MINCOV==TRUE & tabFil$MINEVALUE==TRUE & tabFil$PROPERPAIRS==TRUE & tabFil$MINLEN==TRUE)]),]
# deleteSca <- tabSca[which(tabSca$Scaffold %in% tabFil$ID_SCA[which(tabFil$COVNA==FALSE | tabFil$COVREF==FALSE | tabFil$MINCOV==FALSE | tabFil$MINEVALUE==FALSE | tabFil$PROPERPAIRS==FALSE | tabFil$MINLEN==FALSE)]),]
# 
# 
# ins <- read.table("./scaffolds/K2110059N/ray/ray21/insertHuman/insert.filter.tab", quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE)
# 
# insToDele <- ins[which(ins$ID_SCA %in% deleteSca$Scaffold),]


