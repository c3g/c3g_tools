# echo "R --no-save --args 

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
meanCov <- as.numeric(args[5])
insertSize <- as.numeric(args[6])
strand <- as.numeric(args[7])
filterRead <- as.numeric(args[8])
readLen <- as.numeric(args[9])

#lib

#set path 
folderScaffold <- paste(mainPath,"/scaffolds/", sample,"/ray/ray",kmer,"/",sep="")
folderOut <- paste(folderScaffold,"/insert",type,"/", sep="")
dir.create(file.path(folderOut), showWarnings = FALSE, recursive=TRUE)

# #########################################################################################################
# #########################################################################################################
# FONCTION
# #########################################################################################################
# #########################################################################################################

# #########################################################################################################
# #########################################################################################################
# MAIN
# #########################################################################################################
# #########################################################################################################


file <- paste(folderOut,"cluster.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
allCluster <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "logical", "numeric", "character", "logical"))
allCluster <- allCluster[which(allCluster$ID_INSERT!="-1"),]


file <- paste(folderOut,"insert.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
allInsert <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "character", "character", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))


file <- paste(folderOut,"scaffolds.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
tabSca <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character", "numeric", "numeric"))
colnames(tabSca) <- c("ID_SCA", "NB_READ", "NB_PAIRED_READ", "NB_PROP_PAIRS", "T_COVERED_BASES", "T_COV", "MEAN_COV", "Q75", "Q50", "Q25", "REF_GC", "P_BASES_COV_10X", "P_BASES_COV_100X", "LEN_SCA", "T_MATCH", "MATCH", "NUM_M", "EVALUE")

file <- paste(folderOut,"scaffolds.toDelete.tab",sep="")
if (!file.exists(file) | file.info(file)$size==0){
q()
}
tabScaToDelete <- read.table(file, quote="", comment.char="", stringsAsFactors=FALSE, sep="\t", header=TRUE, colClasses=c("character", "logical", "logical", "logical"))
#deleteSca <- tabScaToDelete[which(tabScaToDelete$COVREF==FALSE | tabScaToDelete$MINCOV==FALSE | tabScaToDelete$PROPERPAIRS==FALSE | tabScaToDelete$MINLEN==FALSE),]
deleteSca <- tabScaToDelete[which(tabScaToDelete$COVREF==FALSE | tabScaToDelete$PROPERPAIRS==FALSE | tabScaToDelete$MINLEN==FALSE),]

# ###################################
# Filter result
# ###################################

#delete bad quality scaffolds 
o <- which(tabSca$ID_SCA %in% deleteSca$ID_SCA)
if (length(o)!=0){
  tabSca <- tabSca[-o,]
}

if (length(tabSca$ID_SCA)==0) {
  q()
}

#delete insert with small cluster
o <- which(abs(allInsert$REF_END-allInsert$REF_START)<=(readLen/2))
if (length(o)!=0){
  allInsert <- allInsert[-o,]
}

if (length(allInsert$ID_SCA)==0) {
  q()
}

#select max insertion by scaffolds
allInsertFilter <- NULL
for (i in 1:length(tabSca$ID_SCA)) {
  nbIns <- (tabSca$MEAN_COV[i]/meanCov)
  if (strand==1) { 
    nbIns <- round(nbIns * 2) 
  }  else { nbIns <- round(nbIns) }
  if (nbIns>0) {
    selectIns <- allInsert[which(allInsert$ID_SCA==tabSca$ID_SCA[i]),]
    if (length(selectIns$ID_INSERT)>0) {
      selectIns <- selectIns[order(selectIns$SCORE,decreasing = TRUE),]
      nb <- nbIns
      if (nb>length(selectIns$ID_INSERT))	{nb <- length(selectIns$ID_INSERT)}
      selectIns <- selectIns[c(1:nb),]

      allInsertFilter <- rbind(allInsertFilter,selectIns)
    }
  }
}

allClusterFilter <- allCluster[which(allCluster$ID_INSERT %in% allInsertFilter$ID_INSERT),]
allInsertFilter <- allInsertFilter[which(allInsertFilter$ID_INSERT %in% allClusterFilter$ID_INSERT),]
tabScaFilter <- tabSca[which(tabSca$ID_SCA %in% allInsertFilter$ID_SCA),]


#Delete insertion with #cluster <= filterRead
del <- which((allInsertFilter$SCLIP1+allInsertFilter$OEA1 + allInsertFilter$SCLIP2+allInsertFilter$OEA2)<=filterRead)
if (length(del)>0) {
  allInsertFilter <- allInsertFilter[-del,]
  allClusterFilter <- allClusterFilter[which(allClusterFilter$ID_INSERT %in% allInsertFilter$ID_INSERT),]
  tabScaFilter <- tabScaFilter[which(tabScaFilter$ID_SCA %in% allInsertFilter$ID_SCA),]
}

#Delete insertion without sclip cluster or OEA cluster
del <- which((allInsertFilter$SCLIP1 + allInsertFilter$SCLIP2)==0 | (allInsertFilter$OEA1+allInsertFilter$OEA2)==0)
if (length(del)>0) {
  allInsertFilter <- allInsertFilter[-del,]
  allClusterFilter <- allClusterFilter[which(allClusterFilter$ID_INSERT %in% allInsertFilter$ID_INSERT),]
  tabScaFilter <- tabScaFilter[which(tabScaFilter$ID_SCA %in% allInsertFilter$ID_SCA),]
}

write.table(tabScaFilter, file=paste(folderOut,"scaffolds.filter.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)
write.table(allInsertFilter, file=paste(folderOut,"insert.filter.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)
write.table(allClusterFilter, file=paste(folderOut,"cluster.filter.tab",sep=""), sep="\t", dec = ".", row.names = FALSE, col.names= TRUE, quote=FALSE)





