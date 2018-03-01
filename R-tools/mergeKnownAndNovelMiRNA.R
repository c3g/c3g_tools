
#function to gather known miRNA and merge them with the list known miRNAs from mirdeep2
#by François Lefebvre, adapted by Eloi Mercier

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
options(stringsAsFactors=F)

# cd $PHOME
# export known_mature="$PHOME/mirbase/mature_target.fa" # need to export so that visible from R
# export known_precursor="$PHOME/mirbase/precursors_target.fa" # need to export so that visible from R
# export odir="$PHOME/mirdeep2"; mkdir -p $odir

# Parse the csv report to extract novel miRNAs
# fn=list.files("mirdeep2",pattern='^result.*csv$',full.names=T)
# odir  = Sys.getenv("odir")
# species = Sys.getenv("MIRBASE_TARGET_SPECIES_CODES")

args=commandArgs(TRUE)
mirdeep_result_file=args[1]
species=args[2]
output_dir=args[3]
# species = Sys.getenv("MIRBASE_TARGET_SPECIES_CODES")

# Parse out the novel table
lines = readLines(fn)
lines = lines[lines!='']
novel = lines[(1+grep("novel miRNAs predicted by miRDeep2",lines)):(grep("mature miRBase miRNAs detected by miRDeep2",lines)-1)]
novel = read.delim(textConnection(novel),check.names=FALSE,stringsAsFactors=F)

# mirbase-style names
novel[["mature.dummy.mirbase.id"]]  = paste0(species,"-miR-"   ,paste0("novel",1:nrow(novel))  )
novel[["precursor.dummy.mirbase.id"]]  = paste0(species,"-mir-", paste0("novel",1:nrow(novel)) )
rownames(novel) = paste(novel[["mature.dummy.mirbase.id"]],"from",novel[["precursor.dummy.mirbase.id"]] ,sep=' ')

# Nice description string (SYMBOL)
novel[["SYMBOL"]] = novel[["mature.dummy.mirbase.id"]]
novel[["SYMBOL"]][ novel[["example miRBase miRNA with the same seed"]] !='-' ] = novel[["example miRBase miRNA with the same seed"]][ novel[["example miRBase miRNA with the same seed"]] !='-' ]

# create a track from the novel (already )
x = strsplit(novel[["precursor coordinate"]],split="(:|\\.\\.)")
x = DataFrame(do.call(rbind,x))
colnames(x) = c("chr","start","end","strand")
x$start = as.numeric(x$start)
x$end = as.numeric(x$end)
x =  makeGRangesFromDataFrame(x)
# export.bed(x,file.path(odir,"novel_miRNAs.bed"))

# Write to disk for future query
write.csv(novel,file=file.path(odir,"novel_miRNAs.csv"),row.names=F)
# save(novel,file=file.path(odir,"novel_miRNAs.RData"))


# Read the known, merge with novel and write out (don't known what to do with star sequences)
all.seqs = sapply(c("mature","precursor"),function(type)
{
	# novels
	n = RNAStringSet(novel[[sprintf("consensus %s sequence",type)]])
	names(n) = novel[[paste0(type,".dummy.mirbase.id")]]

	# known
	k =  readRNAStringSet(Sys.getenv(sprintf("known_%s",type)))

	# merge snd write
	all = c(k,n)
	writeXStringSet( all, file = file.path(odir, sprintf("all_%s.fa",type)    ) )

	return(all)
},simplify=FALSE)
