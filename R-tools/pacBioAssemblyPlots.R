# Function that takes a infile: (filtered_summary.csv) and an outdir where to
# To write plot files.
pacBioPlots <- function(infile, outdir) {
  library(ggplot2)
  
  tData <- read.csv(file=infile, header=T)
  tData <- tData[tData[8] == 1,]
  #readLength <- tData[,4] + 1
  
  outfile1 <- file.path(outdir, "pacBioGraph_readLengthScore.pdf")
  outfile2 <- file.path(outdir, "pacBioGraph_readLengthScore.jpeg")
  outfile3 <- file.path(outdir, "pacBioGraph_histoReadLength.pdf")
  outfile4 <- file.path(outdir, "pacBioGraph_histoReadLength.jpeg")
  
  p <- ggplot(environment = environment(), data=tData, aes(x=tData[,4], y=tData[,5])) + 
    geom_point(shape=20, fill="#0066CC") +   
    xlab("Read length") +
    ylab("Read quality score") +
    theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))
  pdf( file=outfile1 )
  print(p)
  dev.off()  
  
  p <- ggplot(environment = environment(), data=tData, aes(x=tData[,4], y=tData[,5])) + 
    geom_point(shape=20, fill="#0066CC") +   
    xlab("Read length") +
    ylab("Read quality score") +
    theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))
  jpeg( file=outfile2 )
  print(p)
  dev.off()
  
  p <- ggplot(environment = environment(), data=tData, aes(x=readLength)) +
    geom_histogram(binwidth=100, fill="#0066CC") +
    #scale_y_log10() +
    xlab("Read length") +
    ylab("Frequency") +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"))
  pdf( file=outfile3 )
  print(p)
  dev.off()
  
  p <- ggplot(environment = environment(), data=tData, aes(x=readLength)) +
    geom_histogram(binwidth=100, fill="#0066CC") +
    #scale_y_log10() +
    xlab("Read length") +
    ylab("Frequency") +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"))
  jpeg( file=outfile4 )
  print(p)
  dev.off()
} 

usage=function(errM) {
        cat("\nUsage : Rscript pacBioAssemblyPlots.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outdir\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if (ARG[i] == "-i") {
		infile=ARG[i+1]
	} else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}
}

pacBioPlots(infile, outdir)

