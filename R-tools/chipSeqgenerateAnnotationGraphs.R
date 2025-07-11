# Generates graphs based on different annotation statistics
# Maxime Caron - Jan 2012
# Modified by Johanna Sandoval, 2013-11-05
# Modified by Paul Stretenowich - Jul 2021

library(plotrix)
args <- commandArgs(TRUE)
design.file<-args[1]
readset.file<-args[1]
output_dir<-args[2]
designs<-read.table(design.file, header=F, sep="\t", check.names=F)
readsets<-read.table(readset.file, header=T, sep="\t", check.names=F)
samples_hash_table <- new.env(hash=TRUE)

narrow.peaks=FALSE

for (i in 1:nrow(readsets)) {
    sample_name <- readsets$Sample[i]
    mark_name <- readsets$MarkName[i]
    mark_type <- readsets$MarkType[i]
    samples_hash_table[[sample_name]][[mark_name]] <- unique(append(samples_hash_table[[sample_name]][[mark_name]], mark_type))
}

for (sample_name in ls(samples_hash_table)) {
    toPrint<-rbind(c("Sample", "Mark Name", "Number of peaks", "Percent near tss", "Median peak height", "Highest peak", "Lowest peak", "Avg peak width"))
    for (mark_name in ls(samples_hash_table[[sample_name]])) {
        mark_type <- samples_hash_table[[sample_name]][[mark_name]]
        if(mark_type == "N") {
            narrow.peaks=TRUE

            # TSS categories stats
            postscript(paste("graphs/", sample_name, ".", mark_name, "_Misc_Graphs.ps", sep=""), paper="letter", horizontal=T)
            par(mfrow=c(2,2), cex.main=0.6)
            annotationDir=file.path(output_dir, "annotation")
            tss<-paste(annotationDir, "/", sample_name, "/", mark_name, "/", sample_name, ".", mark_name, ".tss.stats.csv", sep="")
            print(paste("NOTICE: processing annotation tss stats files: ", tss, sep=""))

            prefix<-gsub(".tss.stats.csv", "", tss , fixed=T)
            d1<-read.table(tss, header=T, sep=",", check.names=F)
            for(i in 1:nrow(d1)) {
                slices <- c(d1[,1] + d1[,2], d1[,3], d1[,4], d1[,5], d1[,6], d1[,7])
                lbls <- c("gene", names(d1[3:length(d1)]))
                pct <- round(slices/sum(slices)*100)
                lbls <- paste(lbls, "(",pct, sep="") # add percents to labels
                lbls <- paste(lbls,"%)",sep="") # ad % to labels 
                pie(slices,labels=lbls, main=paste("Location analysis of binding sites\nSample: ", sample_name, "\nMark Name: ", mark_name, sep=""))
            }
            # Exon intron stats
            exons<-paste(prefix, ".exon.stats.csv",sep="")
            print(paste("NOTICE: processing annotation exon stats files: ", paste(exons, collapse=","),sep=" "))
            if(file.info(exons)$size == 0 || length(scan(exons)) == 0) {
                d1=data.frame(c(0))
            }
            else {
                d1<-read.table(exons, header=F, sep=",", check.names=F)
            }
            hist(d1[,1], breaks=length(levels(as.factor(d1[,1]))), xlim=c(0,20), main=paste("Distribution of peaks found within exons\nSample: ", sample_name, "\nMark Name: ", mark_name, sep=""), xlab="Exon", ylab="Number of peaks")
            introns=paste(prefix, ".intron.stats.csv",sep="")
            print(paste("NOTICE: processing annotation intron stats files: ", paste(introns, collapse=","),sep=" "))
            if(file.info(introns)$size == 0 || length(scan(introns)) == 0) {
                d2=data.frame(c(0))
            }
            else {
                d2<-read.table(introns, header=F, sep=",", check.names=F)
            }
            hist(d2[,1], breaks=length(levels(as.factor(d2[,1]))), xlim=c(0,20), main=paste("Distribution of peaks found within introns\nSample: ", sample_name, "\nMark Name: ", mark_name, sep=""), xlab="Intron", ylab="Number of peaks")
            # Distance to TSS
            distance=paste(prefix, ".tss.distance.csv",sep="")
            print(paste("NOTICE: processing annotation distance to tss stats files: ", paste(distance, collapse=","),sep=" "))      
            d1<-read.table(distance, header=F, sep=",", check.names=F)
            d1<-subset(d1, d1[,1] > -10000 & d1[,1] < 10000)
            hist(d1[,1], breaks=seq(-10000,10000,1000), main=paste("Distribution of peak distances relative to TSS\nSample: ", sample_name, "\nMark Name: ", mark_name, sep=""), xlab="Distance to TSS (bp)", ylab="Number of peaks")
            dev.off()
            peakDir=file.path(output_dir, "peak_call")
            peakfile<-paste(peakDir, "/", sample_name, "/", mark_name, "/", sample_name, ".", mark_name, "_peaks.narrowPeak", sep="")
            print(paste("NOTICE: processing annotation tss stats files: ", tss, sep=""))
            averagePeakWidth<-NA
            highestPeak<-NA
            lowestPeak<-NA
            medianPeakHeight<-NA
            percentNearTSS<-NA
            prefix=gsub("_peaks.narrowPeak", "", peakfile)
            d1<-read.table(peakfile, header=F, sep="\t", check.names=F)
            nPeaks<-nrow(d1)
            averagePeakWidth<-round(mean(d1[,3]-d1[,2]),0)
            # read summits file
            bedfile<-paste(prefix, "_summits.bed",sep="")
            if (file.exists(bedfile)){
                d1<-read.table(bedfile, header=F, sep="\t", check.names=F)
                highestPeak<-max(d1[,5])
                lowestPeak<-min(d1[,5])
                medianPeakHeight<-median(d1[,5])
            }
            if (file.exists(distance)){
                d1<-read.table(distance, header=F, sep=",", check.names=F)
                d2<-subset(d1, d1[,1]>-1000 & d1[,1]<1000)
                percentNearTSS<-round(nrow(d2)/nrow(d1), 2)*100
            }
            toPrint<-rbind(toPrint, c(sample_name, mark_name, nPeaks, percentNearTSS, medianPeakHeight, highestPeak, lowestPeak, averagePeakWidth))
        }
    }
    if (narrow.peaks){
        write.table(toPrint, paste(output_dir, "/annotation/", sample_name, "/peak_stats.csv", sep=""), row.names=F, col.names=F, quote=F, sep=",")
    }
}
