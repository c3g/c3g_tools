# Needed for miRDeep2 as part of the input files preparation :
# Extract all relevant miRNAs from mirbase and clean the fastas to meet mirdeep2 expectations: A,C,T,U,N, no long names, etc.
# Written by Edouard Henrion - Feb 2017
# Usage : Rscript mirbasePrepareForMirdeep.R -t target_species -o other_species -m mirbase_path -l local_mirbase_path

library(Biostrings)

usage=function(errM) {
    cat("\nUsage : Rscript mirbasePrepareForMirdeep.R [option] <Value>\n")
    cat("       -t        : target_species mirbase code e.g. hsa for home sapien\n")
    cat("       -o        : other_species mirbase code e.g. mmu for mus musculus\n")
    cat("       -m        : path of the main mirbase db\n")
    cat("       -l        : path of the local mirbase in which will be stored the prepared 'target' and 'other' files\n")
    cat("       -h        : this help\n\n")

    stop(errM)
}

ARG = commandArgs(trailingOnly = T)
## default arg values
target_species=""
other_species=""
mirbase_path=""
## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-t") {
        target_species=ARG[i+1]
    } else if (ARG[i] == "-o") {
        other_species=ARG[i+1]
    } else if (ARG[i] == "-m") {
        mirbase_path=ARG[i+1]
    } else if (ARG[i] == "-l") {
        mirbase_local_path=ARG[i+1]
    } else if (ARG[i] == "-h") {
        usage("")
    }
}
## check arg consitency
if (mirbase_path == "") {
    usage("Error : Output directory not specified")
}
if (!(file.exists(mirbase_path))) {
    usage("Error : Design file not found")
}

# Deal with species of interest
species =  unique(unlist(strsplit(target_species, split=',')))
other.species = setdiff(unique(unlist(strsplit(other_species, split=','))), species)

# Read all and clean all upfront
mirbase = sapply(c(paste(mirbase_path,"hairpin.fa",sep="/"), paste(mirbase_path,"mature.fa",sep="/")), function(fn){
    s = readRNAStringSet(fn) # print(uniqueLetters(s))
    s = RNAStringSet(gsub('[^ACGUN]', 'N', s)) # print(uniqueLetters(s))
    # Parse the fasta header and append as metadata to object
    mcols(s) = DataFrame(
        'full.name' = names(s),
        'clean.name' = gsub(' .*', '', names(s)),
        'species' = gsub('-.*', '', names(s))
    )
    # Rename with clean names
    names(s) = mcols(s)$clean.name
    return(s)
}, simplify=FALSE)

# Read mirbase's organism table
orgs = read.delim(paste(mirbase_path,"organisms.txt",sep="/"), check.names=FALSE, stringsAsFactors=F)
all(unlist(sapply(mirbase, function(s)mcols(s)$species)) %in% orgs[["#organism"]]) # quick check if all species there

# Define mirdeep's miRNAs_ref  (few mature rabbit)
writeXStringSet(mirbase[[paste(mirbase_path,"mature.fa",sep="/")]][mcols(mirbase[[paste(mirbase_path,"mature.fa",sep="/")]])$species %in% species], file=paste(mirbase_local_path,'mature_target.fa',sep="/"))

# Define mirdeep's miRNAs_other (HS)
writeXStringSet(mirbase[[paste(mirbase_path,"mature.fa",sep="/")]][mcols(mirbase[[paste(mirbase_path,"mature.fa",sep="/")]])$species %in% other.species], file=paste(mirbase_local_path,'mature_others.fa',sep="/"))

# Define mirdeep's precursors
writeXStringSet(mirbase[[paste(mirbase_path,"hairpin.fa",sep="/")]][mcols(mirbase[[paste(mirbase_path,"hairpin.fa",sep="/")]])$species %in% species], file=paste(mirbase_local_path,'precursors_target.fa',sep="/"))
