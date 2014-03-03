#!/usr/bin/env Rscript

# Function that takes a Tree file and write plots to pdf and jpeg files.
makeplots <- function(filename, outdir, prefix) {

  library(ggplot2)

  outfileJpeg = paste0(outdir, "/", prefix, ".jpeg")
  outfilePdf  = paste0(outdir, "/", prefix, ".pdf")

  data = read.table(filename, sep="\t", header=F, skip=1)
  midPoint = ((data[,3] - data[,2])/2) + data[,2]
  covValue = data[,5]
  refID = data[,1]
  df = data.frame(refID, midPoint, covValue)

  p <- ggplot(environment = environment()) +
  geom_line(data=df, aes(x=midPoint, y=covValue), col="#4169E1") + 
  facet_wrap(~refID, scales="free") +
  xlab("Genome position") +
  ylab("Coverage(X)") +
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
        cat("\nUsage : Rscript pacBioAssemblyPlots.R [option] <Value>\n")
        cat("       -t        : Coverage bed file\n")
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
		data=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}else if (ARG[i] == "-p") {
		prefix=ARG[i+1]
	}
}

makeplots(data, outdir, prefix)

