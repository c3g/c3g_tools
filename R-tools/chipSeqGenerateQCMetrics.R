# process.R

# Modified by Paul Stretenowich - Feb 2021

args <- commandArgs(TRUE)
readset.list <- args[1]
output_dir <- args[2]
readsets <- read.table(readset.list, header=T, sep="\t", check.names=F)
samples_hash_table <- new.env(hash=TRUE)

for (i in 1:nrow(readsets)) {
    sample_name <- readsets$Sample[i]
    mark_name <- readsets$MarkName[i]
    samples_hash_table[[sample_name]] <- unique(append(samples_hash_table[[sample_name]], mark_name))
}

for (sample_name in ls(samples_hash_table)) {
    for (mark_name in samples_hash_table[[sample_name]]) {
        postscript(paste(output_dir, "/graphs/", sample_name, ".", mark_name, "_QC_Metrics.ps", sep=""))
        par(mfrow=c(2,2))

        # Tag count distribution

        countDist<-paste(output_dir, "tags", sample_name, paste(sample_name, mark_name, sep="."), "tagCountDistribution.txt", sep="/")
        d2<-read.table(countDist, header=T, sep="\t")
        d2<-subset(d2, d2[,1]<=7)
        mean<-0
        for(j in 1:nrow(d2)) {
            mean<-mean+(d2[j,2]*d2[j,1])
        }

        barplot(d2[2:10,2], names.arg=d2[2:10,1],ylim=c(0,1),xlab="Number of Tags per Genomic Position", ylab="Fraction of Genomic Positions", main=paste(sample_name, mark_name, "\nTag count: ", round(mean,3), " avg tag per position", sep=" "), col="blue")

        # Tag auto correlation

        countDist<-paste(output_dir, "tags", sample_name, paste(sample_name, mark_name, sep="."), "tagAutocorrelation.txt", sep="/")
        d1<-read.table(countDist, header=T, sep="\t")
        maxH<-max(max(d1[,2]),max(d1[,3]))

        plot(d1[,1], d1[,2], xlab="Distance from reference tag (bp)", ylab="Number of observations", xlim=c(-300,400), ylim=c(0,maxH+500), main=paste(sample_name, mark_name, " \nTag autocorrelation", sep=" "), col="blue", type='l')
        lines(d1[,1], d1[,3], col="red")
        legend("topleft", c("same strand","opposite strand"), col=c("blue", "red"), pch='-', pt.cex=3,bty='n') 

        # Sequence bias

        countDist<-paste(output_dir, "tags", sample_name, paste(sample_name, mark_name, sep="."), "tagFreq.txt", sep="/")
        d1<-read.table(countDist, header=T, sep="\t")
        plot(d1[,1], d1[,2], xlab="Distance from 5' end of tags", ylab="Nucleotide Frequency", col="blue", xlim=c(-50,100), ylim=c(0,0.45), type='l', main=paste(sample_name, mark_name, "\nSequence bias", sep=" "))
        lines(d1[,1],d1[,3], col="red")
        lines(d1[,1],d1[,4], col="yellow")
        lines(d1[,1],d1[,5], col="green")
        abline(v=0)
        for(i in seq(0,0.45,0.05)) {
        lines(c(0,-5),c(i,i))
        }

        legend("bottomright", c("A","C","G","T"), col=c("blue","red","yellow","green"), pch="-",bty='n',pt.cex=3)

        # GC bias

        countGenome<-paste(output_dir, "tags", sample_name, paste(sample_name, mark_name, sep="."), "tagGCcontent.txt", sep="/")
        countIP<-paste(output_dir, "tags", sample_name, paste(sample_name, mark_name, sep="."), "genomeGCcontent.txt", sep="/")
        d1<-read.table(countGenome, header=T, sep="\t")
        d2<-read.table(countIP, header=T, sep="\t")

        maxH<-max(max(d1[,3]),max(d2[,3]))

        plot(d1[,1], d1[,3], col="black", type='l', xlab="GC content of ChIP fragments", ylim=c(0,maxH),ylab="Normalized Fraction of total fragments", main=paste(sample_name, mark_name, "\n GC bias", sep=" "))
        lines(d2[,1],d2[,3], col="red")

        legend("topright", c("Genome", "ChIP sample"), col=c("black","red"), pch="-", bty='n', pt.cex=3)

        dev.off()
    }
}
