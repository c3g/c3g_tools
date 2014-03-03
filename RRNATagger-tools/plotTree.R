#!/usr/bin/env Rscript

# Function that takes a Tree file and write plots to pdf and jpeg files.
makeplots <- function(treeFile, outdir, prefix) {

  library(ape)

  outfileJpeg = paste0(outdir, "/", prefix, ".jpeg")
  outfilePdf  = paste0(outdir, "/", prefix, ".pdf")

  tree = read.tree(treeFile)
  
  jpeg(file=outfileJpeg, height=11, width=8, units="in", res=500)
  plot(tree, cex=1)
  dev.off() 

  pdf(file=outfilePdf, height=16, width=8)
  plot(tree, cex=1)
  dev.off() 
} 

usage=function(errM) {
        cat("\nUsage : Rscript pacBioAssemblyPlots.R [option] <Value>\n")
        cat("       -t        : Tree file\n")
        cat("       -o        : outdir\n")
        cat("       -p        : prefix\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-t") {
		treeFile=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}else if (ARG[i] == "-p") {
		prefix=ARG[i+1]
	}
}

makeplots(treeFile, outdir, prefix)

