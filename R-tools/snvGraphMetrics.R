## merge all sample readStat
## Mathieu bourgey 
## 2013/01/31

##load libraries
library(pheatmap)

args=commandArgs(TRUE)
listFiles=scan(args[1],sep="\n",what='character')
outputBaseName=args[2]
listFileOutBasename=NULL

## what extension should be included in the listFiles
## Summary.table.csv
## Change.rate.by.chromosome.csv
## Changes.by.type.csv
## Effects.by.impact.csv
## Effects.by.functional.class.csv
## Count.by.effects.csv
## Count.by.genomic.region.csv
## Quality.csv
## Coverage.csv
## InDel.lengths.csv
## Base.changes.csv
## TsTv.summary.csv
## TsTv.All.variants.csv
## TsTv.Known.variants.csv
## Allele.frequency.csv
## Allele.frequency.All.variants.csv
## Allele.frequency.Known.variants.csv
## Codon.change.table.csv
## Amino.acid.change.table.csv
## Chromosome.change.table.csv
## changeRate.tsv
##

fileExtensionRetained=c("Summary.table.csv","TsTv.summary.csv","changeRate.tsv","Effects.by.impact.csv","Effects.by.functional.class.csv","Count.by.effects.csv","Count.by.genomic.region.csv","Quality.csv","Coverage.csv","InDel.lengths.csv","Base.changes.csv","TsTv.All.variants.csv","TsTv.Known.variants.csv","Codon.change.table.csv","Amino.acid.change.table.csv","Chromosome.change.table.csv")

## Summary table
##   * Number_of_variants_before_filter
##   * Number_of_variants_filtered_out
##   * Number_of_not_variants
##   * Number_of_variants_processed
##   * Number_of_known_variants
## Ts/Tv summary
##  * should be add to the summary table

sumT=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[1],listFiles)],sep="\n",what='character')),",")
valueSum=c("Number_of_variants_before_filter","Number_of_variants_filtered_out","Number_of_not_variants","Number_of_variants_processed","Number_of_known_variants")
summaryTable=NULL
for (i in 1:length(valueSum)){
	summaryTable=rbind(summaryTable,sumT[[grep(valueSum[i],sumT)]][1:2])
}
summaryTable=gsub("<br>(i.e.non-emptyID)","",summaryTable,fixed=T)
perS=list(c("%",as.character(round((as.numeric(summaryTable[2,2])/as.numeric(summaryTable[1,2]))*100,2))),c("%",as.character(round((as.numeric(summaryTable[3,2])/as.numeric(summaryTable[1,2]))*100,2))),c("%",as.character(round((as.numeric(summaryTable[5,2])/as.numeric(summaryTable[4,2]))*100,2))))
summaryTable=rbind(summaryTable[1:2,],perS[[1]],summaryTable[3,],perS[[2]],summaryTable[4:5,],perS[[3]])
sumTs=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[2],listFiles)],sep="\n",what='character')),",")
for ( i in 1:length(sumTs)) {
	summaryTable=rbind(summaryTable,sumTs[[i]][1:2])
}
colnames(summaryTable)=c("Summary_stats","Value")
write.table(t(summaryTable),paste(outputBaseName,"SummaryTable.tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"SummaryTable",sep="."))

## change rate
##   * graphs
##
chR=read.table(listFiles[grep(fileExtensionRetained[3],listFiles)],header=T,check.names=F,row.names=1)
cn=as.numeric(gsub("chr","",rownames(chR)))
cnpos=1:length(cn)
cnNu=order(cn[cnpos[!(is.na(cn))]])
cnNuV=chR[!(is.na(cn)),]
cnCh=order(rownames(chR)[cnpos[is.na(cn)]])
cnChV=chR[is.na(cn),]
chR.ord=rbind(cnNuV[cnNu,],cnChV[cnCh,])
jpeg(paste(outputBaseName,"changeRate.jpeg",sep="."),800,800)
pheatmap(t(as.matrix(chR.ord)),cluster_cols =F,cluster_rows =F,main="Change rate by sample and by chromosome",fontsize_row=6,fontsize_col=10,display_numbers=T,number_format="%d",fontsize_number=10)
dev.off()
pdf(paste(outputBaseName,"changeRate.pdf",sep="."),title="Change rate by sample and by chromosome",pointsize=5,paper='special')
pheatmap(t(as.matrix(chR.ord)),cluster_cols =F,cluster_rows =F,main="Change rate by sample and by chromosome",fontsize_row=3,fontsize_col=8,display_numbers=T,number_format="%d",fontsize_number=7)
dev.off()
write.table(t(as.matrix(chR.ord)),paste(outputBaseName,"changeRate.tsv",sep="."),quote=F,row.names=T,col.names=T,sep="\t")
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"changeRate",sep="."))

## Effects by impact
## 
effectIlist=strsplit(gsub("%","",gsub(" ","",scan(listFiles[grep(fileExtensionRetained[4],listFiles)],sep="\n",what='character'))),",")
effectITable=NULL
nameEffect=NULL
for (i in 2:length(effectIlist)){
	effectITable=c(effectITable,effectIlist[[i]][2])
	nameEffect=c(nameEffect,effectIlist[[i]][1])
}
effectITable=rbind(nameEffect,effectITable)
write.table(effectITable,paste(outputBaseName,"EffectsImpact.tsv",sep="."),sep="\t",col.names=T,row.names=F,quote=F)
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"EffectsImpact",sep="."))

## Effects by functional class
## 
effectFlist=strsplit(gsub("%","",gsub(" ","",scan(listFiles[grep(fileExtensionRetained[5],listFiles)],sep="\n",what='character'))),",")
effectFTable=NULL
nameEffect=NULL
for (i in 2:length(effectFlist)){
	effectFTable=c(effectFTable,effectFlist[[i]][2])
	nameEffect=c(nameEffect,effectFlist[[i]][1])
}
effectFTable=rbind(nameEffect,effectFTable)
write.table(effectFTable,paste(outputBaseName,"EffectsFunctionalClass.tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"EffectsFunctionalClass",sep="."))

## Count by effects
##  * deleting:
##     CODON_CHANGE_PLUS_CODON_DELETION
##     CODON_CHANGE_PLUS_CODON_INSERTION
##     CODON_DELETION
##     CODON_INSERTION
## 
countByEffect=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[6],listFiles)],sep="\n",what='character')),",")
removeINDEL=c("CODON_CHANGE_PLUS_CODON_DELETION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_INSERTION")
countTable=NULL
countName=NULL
for (i in 2:length(countByEffect)){
	if (!(countByEffect[[i]][1] %in% removeINDEL)){
		countTable=c(countTable,countByEffect[[i]][2])
		countName=c(countName,countByEffect[[i]][1])
	}
}
countTableF=rbind(countName,countTable)
write.table(countTableF,paste(outputBaseName,"CountEffects.tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
jpeg(paste(outputBaseName,"CountEffects.jpeg",sep="."),800,800)
par(las=2)
par(oma=c(10,4,4,1))
barplot(as.numeric(countTable),names.arg=countName,col=rainbow(length(countTable)),main="Total number of variant by effect type")
dev.off()
pdf(paste(outputBaseName,"CountEffects.pdf",sep="."),title="Total number of variant by effect type",paper='special')
par(las=2)
par(oma=c(10,4,4,1))
barplot(as.numeric(countTable),names.arg=countName,col=rainbow(length(countTable)),main="Total number of variant by effect type")
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"CountEffects",sep="."))

## Count by genomic region
## 
countByRegion=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[7],listFiles)],sep="\n",what='character')),",")
regionOrder=c("UPSTREAM","UTR_5_PRIME","SPLICE_SITE_ACCEPTOR","EXON","SPLICE_SITE_DONOR","INTRON","UTR_3_PRIME","DOWNSTREAM","INTERGENIC")
regionName=c("Up","5'","Splicing acceptor","Exon","Splicing donor","intron","3'","Down","Intergenic")
countTable=NULL
countName=NULL
for (i in 2:length(countByRegion)){
		countTable=c(countTable,countByRegion[[i]][2])
		countName=c(countName,countByRegion[[i]][1])
}
orderR=match(regionOrder,countName)
countTableF=rbind(countName[orderR],countTable[orderR])
write.table(countTableF,paste(outputBaseName,"CountRegions.tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
jpeg(paste(outputBaseName,"CountRegions.jpeg",sep="."),800,800)
par(las=2)
par(oma=c(10,4,4,1))
barplot(as.numeric(countTable[orderR]),names.arg=regionName,col=rainbow(length(countTable)),main="Total number of variant by region type")
dev.off()
pdf(paste(outputBaseName,"CountRegions.pdf",sep="."),title="Total number of variant by region type",paper='special')
par(las=2)
par(oma=c(10,4,4,1))
barplot(as.numeric(countTable[orderR]),names.arg=regionName,col=rainbow(length(countTable)),main="Total number of variant by region type")
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"CountRegions",sep="."))

## Quality
##   * change it as a cumulative
quality=t(read.table(listFiles[grep(fileExtensionRetained[8],listFiles)],sep=",",header=T,check.names=F,row.names=1))
write.table(t(quality),paste(outputBaseName,"SNVQuality.tsv",sep="."),sep="\t",col.names=T,row.names=F,quote=F)
jpeg(paste(outputBaseName,"SNVQuality.jpeg",sep="."),800,800)
par(las=2)
par(mar=c(7,7,7,7))
plot(-1000,-10000, xlim=c(0,1000), ylim=c(0,100),main="Variant count by quality",xlab="Variant Quality",ylab="Cumulative variant sum",axes=F)
for (i in seq(0,100,by=10)) {
	abline(h=i,col='grey',lty=3)
}
points(as.numeric(rownames(quality)),cumsum(quality)/sum(quality)*100,type='b',cex=2,pch="*")
points(as.numeric(rownames(quality)),quality/(max(quality)*1.2)*100,type='b',lty=2,col=2,cex=2,pch="*")
axis(2,at=seq(0,100,by=10),labels=seq(0,100,by=10))
axis(1,at=seq(0,1000,by=50),labels=seq(0,1000,by=50))
Map(function(x,y,z) 
  axis(4,at=x,col.axis=y,labels=z,lwd=0,las=1),
  seq(0,100,by=10),
  rep(2,11),
  as.character(round(seq(0,max(quality)*1.2,length=11)))
)
par(las=0)
axis(4,col=2,at=seq(0,100,by=10),labels=F,lty=2,col.ticks=2)
mtext(side = 4, line = 5, "Variant count",col=2)
dev.off()
pdf(paste(outputBaseName,"SNVQuality.pdf",sep="."),title="Variant count by quality",paper='special')
par(las=2)
par(mar=c(7,7,7,7))
plot(-1000,-10000, xlim=c(0,1000), ylim=c(0,100),main="Variant count by quality",xlab="Variant Quality",ylab="Cumulative variant sum",axes=F)
for (i in seq(0,100,by=10)) {
	abline(h=i,col='grey',lty=3)
}
points(as.numeric(rownames(quality)),cumsum(quality)/sum(quality)*100,type='b',cex=2,pch="*")
points(as.numeric(rownames(quality)),quality/max(quality)*1.2*100,type='b',lty=2,col=2,cex=2,pch="*")
axis(2,at=seq(0,100,by=10),labels=seq(0,100,by=10))
axis(1,at=seq(0,1000,by=50),labels=seq(0,1000,by=50))
Map(function(x,y,z) 
  axis(4,at=x,col.axis=y,labels=z,lwd=0,las=1),
  seq(0,100,by=10),
  rep(2,11),
  as.character(round(seq(0,max(quality)*1.2,length=11)))
)
par(las=0)
axis(4,col=2,at=seq(0,100,by=10),labels=F,lty=2,col.ticks=2)
mtext(side = 4, line = 5, "Variant count",col=2)
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"SNVQuality",sep="."))

## Coverage
##   * change it as a cumulative sum and bin the coverage value
coverage=t(read.table(listFiles[grep(fileExtensionRetained[9],listFiles)],sep=",",header=T,check.names=F,row.names=1))
write.table(t(coverage),paste(outputBaseName,"SNVQuality.tsv",sep="."),sep="\t",col.names=T,row.names=F,quote=F)
jpeg(paste(outputBaseName,"SNVCoverage.jpeg",sep="."),800,800)
par(las=2)
par(mar=c(7,7,7,7))
plot(-1000,-10000, xlim=c(0,10000), ylim=c(0,100),main="Variant count by coverage",xlab="Variant Quality",ylab="Cumulative variant sum",axes=F)
for (i in seq(0,100,by=10)) {
	abline(h=i,col='grey',lty=3)
}
points(as.numeric(rownames(coverage)),cumsum(coverage)/sum(coverage)*100,type='b',cex=2,pch="*")
points(as.numeric(rownames(coverage)),(coverage/(max(coverage)*1.2))*100,type='b',lty=2,col=2,cex=2,pch="*")
axis(2,at=seq(0,100,by=10),labels=seq(0,100,by=10))
axis(1,at=seq(0,10000,by=500),labels=seq(0,10000,by=500))
Map(function(x,y,z) 
  axis(4,at=x,col.axis=y,labels=z,lwd=0,las=1),
  seq(0,100,by=10),
  rep(2,11),
  as.character(round(seq(0,max(coverage)*1.2,length=11)))
)
par(las=0)
axis(4,col=2,at=seq(0,100,by=10),labels=F,lty=2,col.ticks=2)
mtext(side = 4, line = 5, "Variant count",col=2)
dev.off()
pdf(paste(outputBaseName,"SNVCoverage.pdf",sep="."),title="Variant count by coverage",paper='special')
par(las=2)
par(mar=c(7,7,7,7))
plot(-1000,-10000, xlim=c(0,10000), ylim=c(0,100),main="Variant count by coverage",xlab="Variant Quality",ylab="Cumulative variant sum",axes=F)
for (i in seq(0,100,by=10)) {
	abline(h=i,col='grey',lty=3)
}
points(as.numeric(rownames(coverage)),cumsum(coverage)/sum(coverage)*100,type='b',cex=2,pch="*")
points(as.numeric(rownames(coverage)),(coverage/(max(coverage)*1.2))*100,type='b',lty=2,col=2,cex=2,pch="*")
axis(2,at=seq(0,100,by=10),labels=seq(0,100,by=10))
axis(1,at=seq(0,10000,by=500),labels=seq(0,10000,by=500))
Map(function(x,y,z) 
  axis(4,at=x,col.axis=y,labels=z,lwd=0,las=1),
  seq(0,100,by=10),
  rep(2,11),
  as.character(round(seq(0,max(coverage)*1.2,length=11)))
)
par(las=0)
axis(4,col=2,at=seq(0,100,by=10),labels=F,lty=2,col.ticks=2)
mtext(side = 4, line = 5, "Variant count",col=2)
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"SNVCoverage",sep="."))

## InDel lengths
##   * barplot
indelL=t(read.table(listFiles[grep(fileExtensionRetained[10],listFiles)],sep=",",header=T,check.names=F,row.names=1))
write.table(t(indelL),paste(outputBaseName,"IndelLength.tsv",sep="."),sep="\t",col.names=T,row.names=F,quote=F)
jpeg(paste(outputBaseName,"IndelLength.jpeg",sep="."),800,800)
par(las=2)
par(mar=c(7,7,7,7))
ba=barplot(t(indelL),axes=F,axisnames =F,main="INDEL count distribution based on their length")
axis(1,at=round(seq(min(ba),max(ba),length=10)),labels=round(seq(min(as.numeric(rownames(indelL))),max(as.numeric(rownames(indelL))),length=10)))
axis(2)
par(las=0)
mtext(side = 1, line = 5, "INDEL length")
mtext(side = 2, line = 5, "Count")
dev.off()
pdf(paste(outputBaseName,"IndelLength.pdf",sep="."),title="INDEL count distribution based on their length",paper='special')
par(las=2)
par(mar=c(7,7,7,7))
ba=barplot(t(indelL),axes=F,axisnames =F,main="INDEL count distribution based on their length")
axis(1,at=round(seq(min(ba),max(ba),length=10)),labels=round(seq(min(as.numeric(rownames(indelL))),max(as.numeric(rownames(indelL))),length=10)))
axis(2)
par(las=0)
mtext(side = 1, line = 5, "INDEL length")
mtext(side = 2, line = 5, "Count")
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"IndelLength",sep="."))

## Base changes
##   * 3d barplot
## 
baseCh=t(read.table(listFiles[grep(fileExtensionRetained[11],listFiles)],sep=",",header=T,check.names=F,row.names=1))
write.table(baseCh,paste(outputBaseName,"BaseChange.tsv",sep="."),sep="\t",col.names=T,row.names=T,quote=F)
colnames(baseCh)=paste("->",colnames(baseCh),sep="")
jpeg(paste(outputBaseName,"BaseChange.jpeg",sep="."),800,800)
barplot(baseCh,beside=T,col=rainbow(4),legend.text=paste("->",colnames(baseCh),sep=""),main="Base Change count",ylim=c(0,round((max(baseCh)*1.2)/100000)*100000),ylab="Count")
dev.off()
pdf(paste(outputBaseName,"BaseChange.pdf",sep="."),title="Base Change count",paper='special')
barplot(baseCh,beside=T,col=rainbow(4),legend.text=paste("->",colnames(baseCh),sep=""),main="Base Change count",ylim=c(0,round((max(baseCh)*1.2)/100000)*100000),ylab="Count")
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"BaseChange",sep="."))

## Ts/Tv : All variants
## Ts/Tv : Known variants
##   * pull the 2 rate together and plot them by samnple as barplot
## 
tsTvA=t(read.table(listFiles[grep(fileExtensionRetained[12],listFiles)],sep=",",header=T,check.names=F,row.names=1))
colnames(tsTvA)=paste("All_Variant",colnames(tsTvA),sep="_")
tsTvK=t(read.table(listFiles[grep(fileExtensionRetained[13],listFiles)],sep=",",header=T,check.names=F,row.names=1))
colnames(tsTvK)=paste("Known_Variant",colnames(tsTvK),sep="_")
tsTv=cbind(tsTvA,tsTvK)
write.table(baseCh,paste(outputBaseName,"TsTv.tsv",sep="."),sep="\t",col.names=T,row.names=T,quote=F)
jpeg(paste(outputBaseName,"TsTv.jpeg",sep="."),dim(tsTv)[1]*12,800)
par(las=2)
par(mar=c(15,3,3,1))
barplot(t(tsTv[,c(3,6)]),beside=T,col=rainbow(2),legend.text=colnames(tsTv)[c(3,6)],ylim=c(0,max(tsTv[,c(3,6)])*1.2),main="Transition-Transversion rate By sample")
dev.off()
pdf(paste(outputBaseName,"tsTV.pdf",sep="."),title="Transition-Transversion rate By sample",width=(dim(tsTv)[1]*12)/72,height=800/72,paper='special')
par(las=2)
par(mar=c(15,3,3,1))
barplot(t(tsTv[,c(3,6)]),beside=T,col=rainbow(2),legend.text=colnames(tsTv)[c(3,6)],ylim=c(0,max(tsTv[,c(3,6)])*1.2),main="Transition-Transversion rate By sample")
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"TsTv",sep="."))

## Codon change table
##   * 3d barplot
## 
codonChLi=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[14],listFiles)],sep="\n",what='character',skip=1)),",")
codonCh=NULL
codonChRN=NULL
for (i in 2:length(codonChLi)) {
	codonCh=rbind(codonCh,as.numeric(codonChLi[[i]][-1]))
	codonChRN=c(codonChRN,codonChLi[[i]][1])
}
colnames(codonCh)=paste(paste("->",codonChLi[[1]][-1],sep="")," ")
rownames(codonCh)=paste(codonChRN," ")
write.table(codonCh,paste(outputBaseName,"codonChange.tsv",sep="."),sep="\t",col.names=T,row.names=T,quote=F)
jpeg(paste(outputBaseName,"codonChange.jpeg",sep="."),800,800)
par(mar=c(15,1,7,15))
pheatmap(codonCh[-1,-1],cluster_cols =F,cluster_rows =F,fontsize_row=11,fontsize_col=11,show_rownames=T,show_colnames=T,main="Codon change table",display_numbers=T,number_format="%d",fontsize_number=10)
dev.off()
pdf(paste(outputBaseName,"codonChange.pdf",sep="."),title="Coddon change Table",paper='special')
par(mar=c(15,1,7,15))
pheatmap(codonCh[-1,-1],cluster_cols =F,cluster_rows =F,fontsize_row=8,fontsize_col=8,show_rownames=T,show_colnames=T,main="Codon change table",display_numbers=T,number_format="%d",fontsize_number=7)
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"codonChange",sep="."))

## Amino acid change table
##   * barplot
## 
aaChLi=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[15],listFiles)],sep="\n",what='character')),",")
aaCh=NULL
aaChRN=NULL
for (i in 2:length(aaChLi)) {
	aaCh=rbind(aaCh,as.numeric(aaChLi[[i]][-1]))
	aaChRN=c(aaChRN,aaChLi[[i]][1])
}
colnames(aaCh)=paste(paste("->",aaChLi[[1]][-1],sep="")," ")
colnames(aaCh)[match("*",colnames(aaCh))]="+"
rownames(aaCh)=paste(aaChRN," ")
rownames(aaCh)[match("*",rownames(aaCh))]="+"
write.table(aaCh,paste(outputBaseName,"AminoAcidChange.tsv",sep="."),sep="\t",col.names=T,row.names=T,quote=F)
jpeg(paste(outputBaseName,"AminoAcidChange.jpeg",sep="."),800,800)
par(mar=c(15,1,7,15))
pheatmap(aaCh[-1*(1:3),-1*(1:3)],cluster_cols =F,cluster_rows =F,fontsize_row=11,fontsize_col=11,show_rownames=T,show_colnames=T,main="Amino acid change table",display_numbers=T,number_format="%d",fontsize_number=10)
dev.off()
pdf(paste(outputBaseName,"AminoAcidChange.pdf",sep="."),title="Amino acid change Table",paper='special')
par(mar=c(15,1,7,15))
pheatmap(aaCh[-1*(1:3),-1*(1:3)],cluster_cols =F,cluster_rows =F,fontsize_row=8,fontsize_col=8,show_rownames=T,show_colnames=T,main="Amino acid change table",display_numbers=T,number_format="%d",fontsize_number=7)
dev.off()
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"AminoAcidChange",sep="."))
 
## Chromosome change table
##   * supllementary figures in a zip
chChgLi=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[16],listFiles)],sep="\n",what='character')),",")
chChg=list()
cp=0
maxB=0
for (i in 2:length(chChgLi)) {
	if (chChgLi[[i]][2] == "Count") {
		cp=cp+1
		chChg[[cp]]=as.numeric(chChgLi[[i]][-1*(1:2)])
		names(chChg[[cp]])=as.numeric(chChgLi[[i-1]][-1*(1:2)])
		names(chChg)[cp]=chChgLi[[i]][1]
		if ((length(chChgLi[[i]])-2) > maxB) {
			maxB=length(chChgLi[[i]])-2
		}
	}
}
cn=as.numeric(gsub("M","25",gsub("Y","24",gsub("X","23",gsub("chr","",names(chChg))))))
cnpos=1:length(cn)
cnNu=order(cn[cnpos[!(is.na(cn))]])
cnNuV=names(chChg)[!(is.na(cn))]
cnCh=order(names(chChg)[cnpos[is.na(cn)]])
cnChV=names(chChg)[is.na(cn)]
chR.ord=c(cnNuV[cnNu],cnChV[cnCh])
pdf(paste(outputBaseName,"chromosomeChange.pdf",sep="."),title="Chromosme change Table",paper='special',onefile=T)
chChgT=matrix(rep("-",maxB*length(chR.ord)),ncol=maxB)
for (i in 1:length(chR.ord)) {
	pos=match(chR.ord[i],names(chChg))
	chChgT[pos,1:length(chChg[[pos]])]=as.numeric(chChg[[pos]])
	plot(chChg[[pos]],type="l",col=2,main=paste("Number of variant change by Mb - chromosom",chR.ord[i]),xlab="position (Mb)",ylab="Change count / Mb")
}
rownames(chChgT)=chR.ord
colnames(chChgT)=paste("Mb_",as.character(1:maxB),sep="")
dev.off()
write.table(aaCh,paste(outputBaseName,"chromosomeChange.tsv",sep="."),sep="\t",col.names=T,row.names=T,quote=F)
zip(zipfile=paste(outputBaseName,"chromosomeChange.zip"),file=c(paste(outputBaseName,"chromosomeChange.pdf",sep="."),paste(outputBaseName,"chromosomeChange.tsv",sep=".")))
listFileOutBasename=c(listFileOutBasename,paste(outputBaseName,"chromosomeChange",sep="."))
write.table(as.data.frame(listFileOutBasename),paste(outputBaseName,"snvGraphMetrics_listFiles.txt",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
