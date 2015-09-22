#!/usr/bin/env Rscript

# Function that takes a Tree file and write plots to pdf and jpeg files.
mergeTables <- function(blastFile, coverageFile, outdir) {
  options(stringsAsFactors = FALSE)
  library(plyr)

  #root         = "/Users/jtrembla/Documents/pacbio_pipeline/"
  #blastFile    = paste0(root, "blast_report.csv.orig")
  #coverageFile = paste0(root, "coverage.bed")
  #outfile      = paste0(root, "blast_report_with_cov.tsv")

  ## Load tables
  # blastFile ="WGA_Sendo_LEV6574/30X/merSize14/blast/blast_report.csv"
  tBlast = readLines(blastFile)
	tBlast = tBlast[!grepl("^#",tBlast)]
  tBlast = read.table( textConnection(tBlast) , sep="\t", comment.char="", header=FALSE, quote="")
  tCoverage = read.table(coverageFile, sep="\t", comment.char="#", header=FALSE, skip=1)

  ## Compute average coverage for each contigs.
  tStats = ddply(tCoverage,~V1,summarise,mean=mean(V5),sd=sd(V5))

  ## Reformat V1
  tStats$V1 = as.character(lapply(strsplit(as.character(tStats$V1), split="\\|"), "[", 1))
  tBlast$V1 = as.character(lapply(strsplit(as.character(tBlast$V1), split="\\|"), "[", 1))

  ## Add coverage value to blast table.
  tBlastMerged = merge(tBlast, tStats, sort=FALSE)
  tBlastMerged$sd = NULL
  tBlastMerged$mean = as.numeric(tBlastMerged$mean)
  tBlastMerged$mean = round(tBlastMerged$mean, digits=0)


  tBlastMerged = data.frame(tBlastMerged$V1, tBlastMerged$V2, tBlastMerged$mean, tBlastMerged$V3, tBlastMerged$V4, tBlastMerged$V5,
                            tBlastMerged$V6, tBlastMerged$V7, tBlastMerged$V8, tBlastMerged$V9, tBlastMerged$V10, tBlastMerged$V11,
                           tBlastMerged$V12, tBlastMerged$V13, tBlastMerged$V14, tBlastMerged$V15, tBlastMerged$V16
  )
  colnames(tBlastMerged) = c( "query id", "subject id", "Coverage(X)", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score", "Reference name", "Kingdom", "Scientific name", "Common name")
  
  write.table(tBlastMerged, file=paste0(outdir, "/blastCov.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

  ## Then remove standard deviation from tStats table and print the table.
  tStats$sd = NULL
  tStats$mean = as.numeric(tStats$mean)
  tStats$mean = round(tStats$mean, digits=0)
  write.table(tStats, file=paste0(outdir, "/contigsCoverage.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

} 

usage=function(errM) {
        cat("\nUsage : Rscript pacBioAssemblyPlots.R [option] <Value>\n")
        cat("       -c        : Coverage bed file\n")
        cat("       -o        : outdir\n")
        cat("       -b        : blast output file\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-b") {
		blastFile=ARG[i+1]
	}else if (ARG[i] == "-c") {
		coverageFile=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}
}

mergeTables(blastFile, coverageFile, outdir)

