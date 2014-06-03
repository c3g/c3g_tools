#!/usr/bin/env Rscript

# Function that takes a Tree file and write plots to pdf and jpeg files.
makeplots <- function(infiles, outdir) {

  #infiles="contigs21.txt,contigs31.txt,contigs41.txt,contigs51.txt,contigs61.txt,contigs81.txt,contigs91.txt"
  
  library(ggplot2)

  prefix = "histogram_contig_length"

  outfileJpeg = paste0(outdir, "/", prefix, ".jpeg")
  outfilePdf  = paste0(outdir, "/", prefix, ".pdf")

  files = unlist(strsplit(infiles, ","))
  df = NULL

  for(i in 1:length(files)){
    print(paste0("Reading file ", files[i]))

    currFile = files[i];
    currTable = read.table(currFile, sep="\t", header=FALSE)
    name = dirname(currFile)
  
    currTable = cbind(currTable, name)
    df = rbind(df, currTable)  
    
    print(paste0("Finished reading file ", files[i]))
  }

  print("Generating histograms...")
  p <- ggplot(data=df, environment = environment(),  aes(x=V2)) +
  geom_histogram(col="#4169E1") + 
  facet_wrap(~name, scales="fixed") +
  xlab("Contig length") +
  ylab("Frequency") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
    axis.text=element_text(size=14, colour="black"),
    axis.title=element_text(size=16),
    panel.grid.major=element_line(colour="black", linetype="dotted"),
    panel.grid.minor=element_blank(),
    panel.background=element_blank()
  )   
  
  jpeg(file=outfileJpeg, height=5, width=8, units="in", res=500)
  plot(p)
  dev.off() 

  pdf(file=outfilePdf, height=5, width=8)
  plot(p)
  dev.off() 
} 

usage=function(errM) {
        cat("\nUsage : Rscript contigsLengthPlot.R [option] <Value>\n")
        cat("       -i        : infiles\n")
        cat("       -o        : outdir\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-i") {
		infiles=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}
}

makeplots(infiles, outdir)

