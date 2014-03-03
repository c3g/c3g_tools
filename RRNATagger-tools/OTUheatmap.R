#!/usr/bin/env Rscript

# Function that takes a OTU table from qiime, out
# To write plot files.
makeplots <- function(OTU_table, outdir, prefix, n) {

  library(pheatmap)

  filename = OTU_table

  data = read.table(filename, sep="\t", header=F, skip=1, comment.char="")
  outfileJpeg = paste0(outdir,'/', prefix, '.jpeg')
  outfilePdf = paste0(outdir,'/', prefix, '.pdf')

  data = as.matrix(data)
  data2=data[2:nrow(data),3:ncol(data)-1]
  colLabels = data[1,]
  rowLabels = data[,1]
  lineages = data[,ncol(data)]
  data3 = data.frame(data2,  stringsAsFactors = FALSE)
  names(data3) = colLabels[3:length(colLabels)-1]
  IDAndLineages = paste(rowLabels, lineages)
  rownames(data3) = IDAndLineages[2:length(IDAndLineages)]
  mat = data.matrix(data3)

  matOrdered = mat[rev(order(rowSums(mat))),]
  
  if(n == 0){
    matOrderedSubset = matOrdered
  }else{
    matOrderedSubset = matOrdered[1:n,]
  }

  pheatmap(matOrderedSubset, file=outfileJpeg, fontsize_row=6, fontsize_col=6,cellwidth=10, cellheight=6)
  pheatmap(matOrderedSubset, file=outfilePdf, fontsize_row=6, fontsize_col=6,cellwidth=10, cellheight=6)  
} 

usage=function(errM) {
        cat("\nUsage : Rscript pacBioAssemblyPlots.R [option] <Value>\n")
        cat("       -t        : OTU table\n")
        cat("       -o        : outdir\n")
        cat("       -p        : prefix\n")
        cat("       -n        : Number of rows to include in heatmap\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-t") {
		OTU_table=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}else if (ARG[i] == "-p") {
		prefix=ARG[i+1]
	}else if (ARG[i] == "-n") {
	  n=ARG[i+1]
	}
}

makeplots(OTU_table, outdir, prefix, n)

