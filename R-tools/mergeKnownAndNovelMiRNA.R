
#function to gather known miRNA and merge them with the list known miRNAs from mirdeep2
#by François Lefebvre, adapted by Eloi Mercier
#change name to parseKnownandNovelMiRNA.R

library(Biostrings)
options(stringsAsFactors=F)

args=commandArgs(TRUE)
mirdeep_result_file=args[1]
species=args[2]
output_dir=args[3]

# Parse out the novel table
lines = readLines(mirdeep_result_file)
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

# Write to disk for future query
write.csv(novel,file=file.path(output_dir,"novel_miRNAs.csv"),row.names=F)

# Read the known, merge with novel and write out (don't known what to do with star sequences)
all.seqs = sapply(c("mature","precursors"),function(type)
{
	# novels
	n = RNAStringSet(novel[[sprintf("consensus %s sequence",type)]])
	names(n) = novel[[paste0(type,".dummy.mirbase.id")]]

	# known
	k =  readRNAStringSet(file.path("mirbase", paste0(type,"_target.fa"))) #check if mirbase exist?

	# merge snd write
	all = c(k,n)
	writeXStringSet( all, file = file.path(output_dir, sprintf("all_%s.fa",type)    ) )

	return(all)
},simplify=FALSE)
