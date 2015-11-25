## merge all sample readStat
## Mathieu bourgey 
## 2013/01/31

args=commandArgs(TRUE)
fileDir=args[1]
outputFile=args[2]

name=NULL
align=NULL
duplicate=NULL
paternFile=c(".sorted.dup.metrics",".insert_size_metrics",".sorted.dup.recal.coverage.tsv")
## get flagstat metrics
listFile=file.path(fileDir,list.files(fileDir,pattern=paste(paternFile[1],"$",sep=""),recursive=T))
sampleNum=length(listFile)
nameSample=strsplit(basename(listFile),".sorted",fixed=T)
for(i in 1:length(listFile)) {
	name=c(name,nameSample[[i]][1])
	stats=scan(file=listFile[i],what="character",sep="\n")
	infoP=grep("^LIBRARY", stats)
	info=strsplit(stats[infoP],"\t")
	pair=strsplit(stats[infoP+1],"\t")
	readExa=as.numeric(pair[[1]][grep("_EXAMINED",info[[1]])])
	align=c(align,readExa[1]+(2*readExa[2]))
	readDup=as.numeric(pair[[1]][grep("_DUPLICATES",info[[1]])])
	duplicate=c(duplicate,readDup[1]+(2*readDup[2])+(2*readDup[3]))
}
pairOrient=rep("NA",sampleNum)
medianInsS=rep("NA",sampleNum)
meanInsS=rep("NA",sampleNum)
averageDev=rep("NA",sampleNum)
standD=rep("NA",sampleNum)
MeanCov=rep("NA",sampleNum)
base10=rep("NA",sampleNum)
base25=rep("NA",sampleNum)
base50=rep("NA",sampleNum)
base75=rep("NA",sampleNum)
base100=rep("NA",sampleNum)
base500=rep("NA",sampleNum)
base1000=rep("NA",sampleNum)
## get insert size metrics
listFile=file.path(fileDir,list.files(fileDir,pattern=paste(paternFile[2],"$",sep=""),recursive=T))
nameSample=strsplit(basename(listFile),".sorted",fixed=T)
if (sampleNum == length(listFile)) {
	sampleNum=length(listFile)
	for(i in 1:length(listFile)) {
		nameLoc=match(nameSample[[i]][1],name)
		stats=scan(file=listFile[i],what="character",sep="\n")
		statsLigne=stats[grep("WIDTH_OF_10_PERCENT",stats,fixed=T)+1]
		statsValue=strsplit(statsLigne,"\t",fixed=T)
		pairOrient[nameLoc]=statsValue[[1]][8]
		medianInsS[nameLoc]=statsValue[[1]][1]
		meanInsS[nameLoc]=statsValue[[1]][5]
		averageDev[nameLoc]=statsValue[[1]][2]
		standD[nameLoc]=statsValue[[1]][6]
	}
} else {
	sampleNumTmp=length(listFile)
	for(i in 1:length(listFile)) {
		if (nameSample[[i]][1] %in% name) {
			nameLoc=match(nameSample[[i]][1],name)
			stats=scan(file=listFile[i],what="character",sep="\n")
			statsLigne=stats[grep("WIDTH_OF_10_PERCENT",stats,fixed=T)+1]
			statsValue=strsplit(statsLigne,"\t",fixed=T)
			pairOrient[nameLoc]=statsValue[[1]][8]
			medianInsS[nameLoc]=statsValue[[1]][1]
			meanInsS[nameLoc]=statsValue[[1]][5]
			averageDev[nameLoc]=statsValue[[1]][2]
			standD[nameLoc]=statsValue[[1]][6]
		} else {
			sampleNum=sampleNum+1
			stats=scan(file=listFile[i],what="character",sep="\n")
			name=c(name,nameSample[[i]][1])
			align=c(align,"NA")
			duplicate=c(duplicate,"NA")
			statsLigne=stats[grep("WIDTH_OF_10_PERCENT",stats,fixed=T)+1]
			statsValue=strsplit(statsLigne,"\t",fixed=T)
			pairOrient=c(pairOrient,statsValue[[1]][8])
			medianInsS=c(medianInsS,statsValue[[1]][1])
			meanInsS=c(meanInsS,statsValue[[1]][5])
			averageDev=c(averageDev,statsValue[[1]][2])
			standD=c(standD,statsValue[[1]][6])
			MeanCov=c(MeanCov,"NA")
			base10=c(base10,"NA")
			base25=c(base25,"NA")
			base50=c(base50,"NA")
			base75=c(base75,"NA")
			base100=c(base100,"NA")
			base500=c(base500,"NA")
                        base1000=c(base1000,"NA")
		}
	}
}


## get coverage metrics
listFile=file.path(fileDir,list.files(fileDir,pattern=paste(paternFile[3],"$",sep=""),recursive=T))
if (length(listFile) > 0) {
	nameSample=strsplit(basename(listFile),".sorted",fixed=T)
	if (sampleNum == length(listFile)) {
		sampleNum=length(listFile)
		for(i in 1:length(listFile)) {
			nameLoc=match(nameSample[[i]][1],name)
			stats=scan(file=listFile[i],what="character",sep="\n", skip=1)	#JT: Skip header, because some samples's name can match header fields...
			statsLigne=stats[grep("Total",stats,fixed=T)]
			statsValue=strsplit(statsLigne,"\t",fixed=T)
			MeanCov[nameLoc]=statsValue[[1]][9]
			base10[nameLoc]=statsValue[[1]][14]
			base25[nameLoc]=statsValue[[1]][15]
			base50[nameLoc]=statsValue[[1]][16]
			base75[nameLoc]=statsValue[[1]][17]
			base100[nameLoc]=statsValue[[1]][18]
			base500[nameLoc]=statsValue[[1]][19]
                        base1000[nameLoc]=statsValue[[1]][20]
		}
	} else {
		sampleNumTmp=length(listFile)
		for(i in 1:length(listFile)) {
			if (nameSample[[i]][1] %in% name) {
				nameLoc=match(nameSample[[i]][1],name)
				stats=scan(file=listFile[i],what="character",sep="\n")
				statsLigne=stats[grep(nameSample[[i]][1],stats,fixed=T)]
				statsValue=strsplit(statsLigne,"\t",fixed=T)
                                MeanCov[nameLoc]=statsValue[[1]][9]
                                base10[nameLoc]=statsValue[[1]][14]
                                base25[nameLoc]=statsValue[[1]][15]
                                base50[nameLoc]=statsValue[[1]][16]
                                base75[nameLoc]=statsValue[[1]][17]
                                base100[nameLoc]=statsValue[[1]][18]
                                base500[nameLoc]=statsValue[[1]][19]
                                base1000[nameLoc]=statsValue[[1]][20]
			} else {
				sampleNum=sampleNum+1
				stats=scan(file=listFile[i],what="character",sep="\n")
				name=c(name,nameSample[[i]][1])
				align=c(align,"NA")
				duplicate=c(duplicate,"NA")
				pairOrient=c(pairOrient,"NA")
				medianInsS=c(medianInsS,"NA")
				meanInsS=c(meanInsS,"NA")
				averageDev=c(averageDev,"NA")
				standD=c(standD,"NA")
                                MeanCov[nameLoc]=statsValue[[1]][9]
                                base10[nameLoc]=statsValue[[1]][14]
                                base25[nameLoc]=statsValue[[1]][15]
                                base50[nameLoc]=statsValue[[1]][16]
                                base75[nameLoc]=statsValue[[1]][17]
                                base100[nameLoc]=statsValue[[1]][18]
                                base500[nameLoc]=statsValue[[1]][19]
                                base1000[nameLoc]=statsValue[[1]][20]
			}
		}
	}
} 
finalTable=cbind(name,align,align - duplicate,duplicate,duplicate / align * 100,pairOrient,medianInsS,meanInsS,averageDev,standD,wgMeanCov,wgbase10,wgbase25,wgbase50,wgbase75,wgbase100,wgbase500,ccdsMeanCov,ccdsbase10,ccdsbase25,ccdsbase50,ccdsbase75,ccdsbase100,ccdsbase500)
colnames(finalTable)=c("Sample","Mapped Reads","Not Duplicate Reads","Duplicate Reads","Duplicate %","Pair Orientation","Median Insert Size","Mean Insert Size","Average Deviation","Standard Deviation","Mean Coverage","%_bases_above_10","%_bases_above_25","%_bases_above_50","%_bases_above_75","%_bases_above_100","%_bases_above_500")
write.table(finalTable,file=outputFile,sep="\t",row.names=F,col.names=T,quote=F)
