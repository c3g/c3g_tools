# Needed for miRDeep2 as aprt of preparation of input files :
# Take a list of arf files and the correponding list of fasta files and filter for mimimim read lengh 17
# Written by Edouard Henrion - March 2017
# Usage : Rscript filterForMirdeep.R -arf arf_files -fasta fasta_files

library(Biostrings)

usage=function(errM) {
    cat("\nUsage : Rscript mirbasePrepareForMirdeep.R [option] <Value>\n")
    cat("       -arf      : list of arf_files\n")
    cat("       -fasta    : list of fasta_files\n")
    cat("       -h        : this help\n\n")

    stop(errM)
}

ARG = commandArgs(trailingOnly = T)
## default arg values
arf_files=""
fasta_files=""
mirbase_path=""
## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-arf") {
        arf_files=ARG[i+1]
    } else if (ARG[i] == "-fasta") {
        fasta_files=ARG[i+1]
    } else if (ARG[i] == "-h") {
        usage("")
    }
}
## check arg consitency
if (arf_files == "") {
    usage("Error : Input arf file list is empty or was not specified")
}
if (fasta_files == "") {
    usage("Error : Input fasta file list is empty or was not specified")
}

arf = read.delim(arf_files, header=F, stringsAsFactors=F, colClasses='character')
arf = arf[ as.numeric(arf[,2]) >= 17,]
write.table(arf, file=arf_files, sep='\t', quote=F, col.names=F, row.names=F)
rm(arf);gc()
fa = readDNAStringSet(fasta_files)
fa = fa[width(fa) >= 17]
writeXStringSet(fa,file=fasta_files)
