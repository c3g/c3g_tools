#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)
#print(args)
suppressPackageStartupMessages(library(gqSeqUtils))
x = calculateGeneLengthsFromGtf(gtf=args[1], as.data.frame=TRUE, output.file=args[2])

#R --no-restore --no-save<<'EOF'
#library(gqSeqUtils)
#x=calculateGeneLengthsFromGtf("annotations/transcripts_ensembl.gtf",feature.types=c("exon","CDS"), as.data.frame=TRUE, output.file="temp.tsv")
#EOF
 
