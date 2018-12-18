## FusionMetaCaller 
## Require multiple fusion detection programs preprocessed - fusion<tab>split+discordant reads
## Robert Eveleigh - robert.eveleigh@mcgill.ca
## 2018/06/18
## usage:
##   Rscript RunFusionMetaCaller.R <INPUT_FOLDER> <OUTPUT_BASE_NAME>

##get args
args=commandArgs(TRUE)
inputFolder=args[1]
sampleName=args[2]

outputFile1 = paste(sampleName, "metacaller.1.tsv", sep=".")
outputFile2 = paste(sampleName, "metacaller.2.tsv", sep=".")
outputFile3 = paste(sampleName, "fusions.tsv", sep=".")

library(FusionMetaCaller)
library(magrittr)

fns<-list.files(path=inputFolder, pattern="*.tsv", full.names=TRUE )

fns

l = sapply(fns,function(fn){
    x = read.table(fn, sep="\t")
    split<-strsplit(fn,".", fixed=TRUE)[[1]]
    colnames(x)<-c("fusion",split[2])
    x
},simplify=FALSE)

res <- Reduce(function(x,y) {merge(x,y,all=T, by="fusion")}, l)
res[is.na(res)]<-0

rownames(res)=res[,1] %>% as.character %>%make.unique

#unique(res)
#rownames(res)=res[,1] 

res = res[,-1]
as.matrix(res)

fmc1<-FusionMetaCaller(res,1)
fmc2<-FusionMetaCaller(res)

fusions<-rownames(fmc2$sortMatrix)

write.table(fmc1, file=outputFile1, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(fmc2, file=outputFile2, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#write.table(as.data.frame(fusions), file=outputFile3, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
