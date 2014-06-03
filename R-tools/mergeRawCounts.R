#!/usr/bin/env Rscript

# Function that takes many raw count matrix table typically generated
# by the MUGQIC RNASeq pipeline and merge them.
# Author: Julien Tremblay - julien.tremblay@mail.mcgill.ca
mergeTables <- function(infiles, outfile) {

  files = unlist(strsplit(infiles, ","))

  df = NULL
  
  tData = read.table(files[1], sep="\t", header=TRUE)

  symbol = as.character(tData[,1])
  geneName = as.character(tData[,2])

  df = cbind(df, symbol)
  df = cbind(df, geneName)

  for(i in 1:length(files)){
    print(paste0("Reading file ", files[i]))
    currFile = files[i];
    currTable = read.table(currFile, sep="\t", header=TRUE)
  
    df = cbind(df, currTable[,3:ncol(currTable)])
  }
	
  write.table(df, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)

	print("Done mergingTables...")

} 

usage=function(errM) {
        cat("\nUsage : Rscript mergeRawCounts.R [option] <Value>\n")
        cat("       -i        : List of files separated by a comma\n")
        cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if (ARG[i] == "-i") {
		infiles=ARG[i+1]
	} else if (ARG[i] == "-o") {
		outfile=ARG[i+1]
	}
}

mergeTables(infiles, outfile)

