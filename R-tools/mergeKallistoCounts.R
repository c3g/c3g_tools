#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

abundance_transcripts_files=args[1]
output_dir=args[2]
data_type=args[3]
#abundance_transcripts_files="kallisto/KW390ASoy200/abundance_genes.tsv,kallisto/KW425Soy10/abundance_genes.tsv"
# output_dir="kallisto"

files = unlist(strsplit(abundance_transcripts_files, ","))
sample_names=sapply(files, function(x) {rev(strsplit(dirname(x),"/")[[1]])[1]})
output_file_name=paste0("all_samples.abundance_" ,data_type, ".csv")

if (data_type=="trasncripts"){
	key_id="target_id"
	count_id="est_counts"
} else if (data_type=="genes") {
	key_id="gene_id"
	count_id="counts"
} else {
	stop("Type must be set to transcrips or genes!")
}
#read first file
f1_table=read.delim(files[1], header=T)
dt_gene_counts=data.frame(matrix(nrow=nrow(f1_table), ncol=length(files)+1))
colnames(dt_gene_counts)=c("geneID",sample_names)
dt_gene_counts[,"geneID"]=f1_table[,key_id]
rownames(dt_gene_counts)=f1_table[,key_id]

for (i in seq(length(files))){
	f_table=read.delim(files[i], header=T, row.names=1)
	dt_gene_counts[,i+1]=f_table[rownames(dt_gene_counts),count_id]
}

write.table(x=dt_gene_counts, file=file.path(output_dir,output_file_name), col.names=T, row.names=F, sep="\t", quote=F)
