library(methylKit)
library(graphics)
library(GenomicRanges)
library(IRanges)

# Usage
usage = function(errM) {
    cat("\nUsage : Rscript deseq.R [option] <Value>\n")
    cat("       -design         : design file\n")
    cat("       -outdir         : output directory\n")
    cat("       -build          : genome build\n")
    cat("       -suff           : input files suffix (inputs should be name as <SAMPLE><suffix>, SAMPLE being the same as in the design file)\n")
    cat("       -mread          : minimm reads\n")
    cat("       -tile           : tile size\n")
    cat("       -step           : step size\n")
    cat("       -mingc          : minimum CpG sites\n")
    cat("       -minmethdiff    : minimum methylation difference\n")
    cat("       -mergestr       : merge strand (True or False)\n")
    cat("       -minsamplegrp   : minimum sample group size\n")
    cat("       -adjust         : adjust option (default 'SLIM')\n")
    cat("       -qvalue         : qvalue threshold (default 0.01)\n")
    cat("       -overdisp       : overdispersion option : 'MN', 'shrinkMN', 'none' (default 'none')\n")
    cat("       -h              : this help\n\n")

    stop(errM)
}

ARG <- commandArgs(trailingOnly=T)

## default arg values
designFile = ""
output_folder = ""
genomeVersion = ""
minReads = 10
tileSize = 500
stepSize = 500
minCGs = 3
minMethDiff = 10
mergeStrand = FALSE
minSamplePerGroup = NULL
qvalue_threshold = 0.01
overdispersion_option = "none"
adjust_option = "SLIM"

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-design") {
        designFile = ARG[i+1]
    } else if (ARG[i] == "-outdir") {
        output_folder = ARG[i+1]
    } else if (ARG[i] == "-build") {
        genomeVersion = ARG[i+1]
    } else if (ARG[i] == "-suff") {
        inputFileSuffix = ARG[i+1]
    } else if (ARG[i] == "-mread") {
        minReads = as.numeric(ARG[i+1])
    } else if (ARG[i] == "-tile") {
        tileSize = as.numeric(ARG[i+1])
    } else if (ARG[i] == "-step") {
        stepSize = as.numeric(ARG[i+1])
    } else if (ARG[i] == "-mcpg") {
        minCGs = as.numeric(ARG[i+1])
    } else if (ARG[i] == "-mmethdiff") {
        minMethDiff = as.numeric(ARG[i+1])
    } else if (ARG[i] == "-mergestr") {
        mergeStrand = as.logical(ARG[i+1])
    } else if (ARG[i] == "-msamplegrp") {
        minSamplePerGroup = as.integer(ARG[i+1])
    } else if (ARG[i] == "-adjust") {
        adjust_option = ARG[i+1]
    } else if (ARG[i] == "-qvalue") {
        qvalue_threshold = as.numeric(ARG[i+1])
    } else if (ARG[i] == "-overdisp") {
        overdispersion_option = ARG[i+1]
    } else if (ARG[i] == "-h") {
        usage("")
    }
}

# if (minSamplePerGroup <= 3) {
#     overdispersion_option <- "none";
# }

## check arg consitency
if (!(file.exists(designFile))) {
    usage("Error : Design file not found !!!")
}
if (output_folder == "") {
    usage("Error : Output directory not specified !!!")
}
if (genomeVersion == "") {
    usage("Error : Genome Build not specified !!!")
}

print(minSamplePerGroup)

design <- read.table(designFile, header=T, stringsAsFactors=F)
##0: skipped sample; 1: contol sample; 2: case sample.

# Iterate per design
for (i in 2:ncol(design)) {
    designName <- colnames(design)[i]
    dir <- paste("mkdir -p ", output_folder, "/graphs/", designName, sep = "")
    system(dir)
    dir <- paste("mkdir -p ", output_folder, "/dmr_results/", designName, sep="")
    system(dir)
    dir <- paste("mkdir -p ", output_folder, "/additional_files/", designName, sep="")
    system(dir)
    dir <- paste("mkdir -p ", output_folder, "/Rdata_files/", designName, sep="")
    system(dir)

    condition1SampleNames <- design[design[,i] == 1,1] #unlist(strsplit(design[i,2],","))
    condition2SampleNames <- design[design[,i] == 2,1] #unlist(strsplit(design[i,3],","))
    condition1Samples <- paste("methylkit/inputs/", condition1SampleNames, inputFileSuffix, sep="")
    condition2Samples <- paste("methylkit/inputs/", condition2SampleNames, inputFileSuffix, sep="")
    nmbSamples <- length(condition1Samples) + length(condition2Samples)
    sampleNames <- c(condition1SampleNames, condition2SampleNames)
    file.list <- as.list(c(condition1Samples, condition2Samples))
    conditions <- c(rep(1, length(condition1Samples)), rep(0, length(condition2Samples)))
    print(paste(sampleNames, collapse=" "))
    print(paste(conditions, collapse=" "))

    myobj <- methRead(file.list, sample.id = as.list(sampleNames), assembly = genomeVersion, treatment = conditions, context = "CpG")
    print("Finishing reading sample methylation profile files!");
    filtered.myobj <- filterByCoverage(myobj, lo.count = minReads, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
    filtered.myobj <- normalizeCoverage(filtered.myobj, "median")
    meth <- unite(filtered.myobj, destrand = mergeStrand,min.per.group = minSamplePerGroup)
    tiles <- tileMethylCounts(myobj, win.size = tileSize, step.size = stepSize, cov.bases=minCGs)
    filtered.tiles <- filterByCoverage(tiles, lo.count = minReads, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
    filtered.tiles <- normalizeCoverage(filtered.tiles, "median")
    meth.tiles <- unite(filtered.tiles, destrand = mergeStrand,min.per.group = minSamplePerGroup)

    # Write sites of merging samples
    print("Write merged methylation profile per base");
    write.meth <- getData(meth)
    for(j in 1:nmbSamples) {
            colnames(write.meth)[5+3*(j-1)] <- paste("coverage_", sampleNames[j], sep="")
            colnames(write.meth)[6+3*(j-1)] <- paste("nmb_Cs_", sampleNames[j], sep="")
             colnames(write.meth)[7+3*(j-1)] <- paste("nmb_Ts_", sampleNames[j], sep="")
            write.meth <- cbind(write.meth, as.data.frame(write.meth[,6+3*(j-1)] / (write.meth[,7+3*(j-1)]+write.meth[,6+3*(j-1)])) * 100)
            write.meth[,ncol(write.meth)] <- round(write.meth[,ncol(write.meth)],2)
            colnames(write.meth)[ncol(write.meth)] <- paste("%_meth_",sampleNames[j],sep="")
    }
    # Write merge sites for all samples
    outDir <- paste(output_folder, "/additional_files/", designName, sep="")
    system(paste("mkdir -p ", outDir))
    nCpGs <- nrow(write.meth[,c(1,2,3,seq(5,ncol(write.meth),1))])
    print(nCpGs);
    write.table(write.meth[,c(1,2,3,seq(5,ncol(write.meth),1))], paste(outDir,"/methsites.merged.normalized.csv",sep=""), quote=F, row.names=F, col.names=T,sep="\t")
    print("Finish writing merged methylation profile per base!");

    # Write tiles of merging samples
    print("Write merged methylation profile per tile");
    write.tiles <- getData(meth.tiles)
    for(j in 1:nmbSamples) {
        colnames(write.tiles)[5+3*(j-1)] <- paste("coverage_",sampleNames[j],sep="")
        colnames(write.tiles)[6+3*(j-1)] <- paste("nmb_Cs_",sampleNames[j],sep="")
        colnames(write.tiles)[7+3*(j-1)] <- paste("nmb_Ts_",sampleNames[j],sep="")
        write.tiles <- cbind(write.tiles, as.data.frame(write.tiles[,6+3*(j-1)] / (write.tiles[,7+3*(j-1)]+write.tiles[,6+3*(j-1)])) * 100)
        write.tiles[,ncol(write.tiles)] <- round(write.tiles[,ncol(write.tiles)],2)
        colnames(write.tiles)[ncol(write.tiles)] <- paste("%_meth_",sampleNames[j],sep="")
    }
    # Write merge tiles for all samples
    write.table(write.tiles[,c(1,2,3,seq(5,ncol(write.tiles),1))], paste(outDir,"/methtiles.merged.normalized.csv",sep=""), quote=F, row.names=F, col.names=T,sep="\t")

    tmp <- read.table(paste(output_folder, "/additional_files/",designName,"/methtiles.merged.normalized.csv",sep=""), header=T, sep="\t", stringsAsFactors=F, check.names=F)
    a <- c(paste(tmp[,1],".",tmp[,2],".",tmp[,3],sep=""))
    b <- cbind(tmp[,1:3], a, a, "+", tmp[,4:ncol(tmp)])
    colnames(b)[4:6] <- c("peak","peak","strand")
    b[,2] <- b[,2]-1
    write.table(b, paste(output_folder, "/additional_files/",designName,"/methtiles.merged.normalized.more_columns.bed",sep=""), quote=F, row.names=F, col.names=T, sep="\t")
    print("Finish writing merged methylation profile per tile!");

    # Generate graphs
    graphOut <- paste(output_folder, "/graphs/", designName, "/percent_methylation_distribution.ps", sep="")
    postscript(graphOut)
    for(j in 1:nmbSamples) {
        getMethylationStats(myobj[[j]], plot = T, both.strands = F)
    }
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/coverage_per_base.ps", sep="")
    postscript(graphOut)
    for(j in 1:nmbSamples) {
        getCoverageStats(myobj[[j]], plot = T, both.strands = F)
    }
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/correlation_per_base.ps", sep="")
    postscript(graphOut)
    getCorrelation(meth, plot = T)
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/correlation_per_tile.ps", sep="")
    postscript(graphOut)
    getCorrelation(meth.tiles, plot = T)
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/cluster_per_base.ps", sep="")
    postscript(graphOut)
    clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/cluster_per_tile.ps", sep="")
    postscript(graphOut)
    clusterSamples(meth.tiles, dist = "correlation", method = "ward", plot = TRUE)
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/pca_variance_per_base.ps", sep="")
    postscript(graphOut)
    PCASamples(meth,screeplot=TRUE)
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/pca_variance_per_tile.ps", sep="")
    postscript(graphOut)
    PCASamples(meth.tiles,screeplot=TRUE)
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/pca_per_base.ps", sep="")
    postscript(graphOut)
    PCASamples(meth)
    dev.off()

    graphOut <- paste(output_folder, "/graphs/", designName, "/pca_per_tile.ps", sep="")
    postscript(graphOut)
    PCASamples(meth.tiles)
    dev.off()

    # Generate methyl % graph per tile per sample
    methyl <- write.tiles[,(ncol(write.tiles)-nmbSamples+1):ncol(write.tiles)]
    m <- matrix(nrow=5,ncol=nmbSamples)
    nmbRows <- nrow(methyl)
    for(j in 1:nmbSamples) {
        m[5,j] <- nrow(subset(methyl,methyl[,j]>=80)) / nmbRows
        m[4,j] <- nrow(subset(methyl,methyl[,j]>60 & methyl[,j] < 80)) / nmbRows
        m[3,j] <- nrow(subset(methyl,methyl[,j]>40 & methyl[,j] < 60)) / nmbRows
        m[2,j] <- nrow(subset(methyl,methyl[,j]>20 & methyl[,j] < 40)) / nmbRows
        m[1,j] <- nrow(subset(methyl,methyl[,j]<=20)) / nmbRows
    }
    graphOut <- paste(output_folder, "/graphs/", designName, "/fraction_tiles_methyl.ps", sep="")
    postscript(graphOut)
    barplot(m,col=c(6,5,4,3,2), cex.axis=0.8, names.arg=sampleNames, main="Fraction of 100bp tiles for each 20% methylation ratio bracket",xlab="Samples", ylab="Fraction of tiles")
    dev.off()

    # Changing tiles between samples (methyl diff > 0.2 and significant t-test)
    print("Start to calcualte differentially methylated sites!")
        meth <- select(meth,1:1000) ## for testing purpose
        meth.tiles <- select(meth.tiles,1:1000)

    myDiff <- calculateDiffMeth(meth, num.cores=2, overdispersion=overdispersion_option, adjust=adjust_option) # overdispersion="shrinkMN", adjust="qvalue"
    myDiff.tiles <- calculateDiffMeth(meth.tiles, num.cores=2, overdispersion=overdispersion_option, adjust=adjust_option)

    #save(meth,myDiff,filtered.myobj,file <- paste(output_folder, "/Rdata_files/", designName, "/methylKit.DMC.rda",sep=""))
    #save(meth.tiles,myDiff.tiles,filtered.tiles,file <- paste(output_folder, "/Rdata_files/", designName, "/methylKit.tiles.DMR.rda",sep=""))

    ## write out the differential DMCs/DMRs
    myDiff.hyper <- getMethylDiff(myDiff, difference = minMethDiff, qvalue = qvalue_threshold, type = "hyper")
    myDiff.tiles.hyper <- getMethylDiff(myDiff.tiles, difference = minMethDiff, qvalue = qvalue_threshold, type = "hyper")
    myDiff.hypo <- getMethylDiff(myDiff, difference = minMethDiff, qvalue = qvalue_threshold, type = "hypo")
    myDiff.tiles.hypo <- getMethylDiff(myDiff.tiles, difference =  minMethDiff, qvalue = qvalue_threshold, type = "hypo")

    print("write dmc hypo base results")
    if(nrow(myDiff.hypo)>0) {
        w <- cbind( getData( myDiff.hypo )[, c(1,2,3)], paste( getData(myDiff.hypo )[,1], ".", getData( myDiff.hypo )[,2], ".", getData( myDiff.hypo )[,3], sep=""), getData( myDiff.hypo )[,c(7,4,5,6)])
        colnames(w)[4] <- "dmr.id"
        write.table(w, paste(output_folder, "/dmr_results/", designName, "/hypo.perbase.", minMethDiff, ".qv", qvalue_threshold, ".percent", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
        w[,2] <- w[,2] - 1
        write.table(w, paste(output_folder, "/dmr_results/", designName, "/hypo.perbase.", minMethDiff, ".qv", qvalue_threshold, ".percent.bed", sep=""), quote=F, row.names=F, col.names=F,sep="\t")

    }
    print("write dmr hypo tiles results")
    if(nrow(myDiff.tiles.hypo)>0) {
        w <- cbind(getData( myDiff.tiles.hypo )[,c(1,2,3)], paste(getData(myDiff.tiles.hypo)[,1], ".", getData(myDiff.tiles.hypo)[,2], ".", getData(myDiff.tiles.hypo)[,3],sep=""), getData(myDiff.tiles.hypo)[,c(7,4,5,6)])
        colnames(w)[4] <- "dmr.id"
        write.table(w, paste(output_folder, "/dmr_results/", designName, "/hypo.pertile.", minMethDiff, ".qv", qvalue_threshold, ".percent", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
        w[,2] <- w[,2]-1
        write.table(w, paste(output_folder, "/dmr_results/", designName, "/hypo.pertile.", minMethDiff, ".qv", qvalue_threshold, ".percent.bed", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
    }
    print("write dmc hyper base results")
    if(nrow(myDiff.hyper)>0) {
        w <- cbind( getData( myDiff.hyper )[,c(1,2,3)], paste( getData( myDiff.hyper )[,1], ".", getData( myDiff.hyper )[,2], ".", getData( myDiff.hyper )[,3], sep=""), getData( myDiff.hyper )[, c(7,4,5,6)])
        colnames(w)[4] <- "dmr.id"
        write.table(w, paste(output_folder, "/dmr_results/",designName,"/hyper.perbase.", minMethDiff, ".qv",qvalue_threshold, ".percent",sep=""), quote=F, row.names=F, col.names=F,sep="\t")
        w[,2] <- w[,2]-1
        write.table(w, paste(output_folder, "/dmr_results/", designName, "/hyper.perbase.", minMethDiff, ".qv",qvalue_threshold, ".percent.bed", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
    }
    print("write dmr hyper tiles results")
    if(nrow(myDiff.tiles.hyper)>0) {
        w <- cbind( getData( myDiff.tiles.hyper )[, c(1,2,3)], paste(getData( myDiff.tiles.hyper )[,1], ".", getData( myDiff.tiles.hyper )[,2], ".", getData( myDiff.tiles.hyper )[, 3], sep=""), getData( myDiff.tiles.hyper )[, c(7,4,5,6)])
        colnames(w)[4] <- "dmr.id"
        write.table(w, paste(output_folder, "/dmr_results/", designName, "/hyper.pertile.", minMethDiff, ".qv",qvalue_threshold, ".percent", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
        w[,2] <- w[,2]-1
        write.table(w, paste(output_folder, "/dmr_results/", designName, "/hyper.pertile.", minMethDiff, ".qv",qvalue_threshold, ".percent.bed",sep=""), quote=F, row.names=F, col.names=F,sep="\t")
    }
    
    #CpGsoutputfile <- paste(output_folder, "/Rdata_files/", designName, "/perbase.testingresults.txt",sep="")
    CpGsoutputgzfile <- gzfile(paste(output_folder, "/Rdata_files/", designName, "/perbase.testingresults.txt.gz",sep=""))
    myDiffData <- getData(myDiff)
    myDiffData[,3] <- myDiffData[,2]+1;
    myDiffData[,4] <- NULL;
    myDiffData[,6] <- round( myDiffData[,6],3)
    myDiffData[,4] <- signif(myDiffData[,4],4)
    myDiffData[,5] <- signif(myDiffData[,5],4)
    colnames(myDiffData) <- c("#chr","start","end","pvalue","qvalue","meth.diff");
    myDiffData[,2] <- as.integer(myDiffData[,2]);
    myDiffData[,3] <- as.integer(myDiffData[,3]);
#    write.table(myDiffData[,c(1,2,3,6,4,5)],CpGsoutputfile,quote=F,row.names=F, col.names = T, sep="\t")
    write.table(myDiffData[,c(1,2,3,6,4,5)],CpGsoutputgzfile,quote=F,row.names=F, col.names = T, sep="\t")

    #CpGtilesoutputfile <- paste(output_folder, "/Rdata_files/", designName, "/pertile.testingresults.txt",sep="")
    CpGtilesoutputgzfile <- gzfile(paste(output_folder, "/Rdata_files/", designName, "/pertile.testingresults.txt.gz",sep=""))
    myDiffData.tiles <- getData(myDiff.tiles)
    myDiffData.tiles[,4] <- NULL;
    colnames(myDiffData.tiles) <- c("#chr","start","end","pvalue","qvalue","meth.diff");
    myDiffData.tiles[,2] <- as.integer(myDiffData.tiles[,2]);
    myDiffData.tiles[,3] <- as.integer(myDiffData.tiles[,3]);
    myDiffData.tiles[,6] <- round( myDiffData.tiles[,6],3)
    myDiffData.tiles[,4] <- signif(myDiffData.tiles[,4],4)
    myDiffData.tiles[,5] <- signif(myDiffData.tiles[,5],4)
    write.table(myDiffData.tiles[,c(1,2,3,6,4,5)],CpGtilesoutputgzfile,quote=F,row.names=F, col.names = T, sep="\t")
    print ("END of Differential Methylation Analysis!")
}
