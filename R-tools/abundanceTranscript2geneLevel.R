#!/usr/bin/env Rscript

#Transforms a abundance file from Kallisto at the transcripts level to gene level
args = commandArgs(trailingOnly=TRUE)

abundance_transcripts_file=args[1]
tx2gene_file=args[2]

require(tximport)

index_path=file.path(dirname(gtf_file),"cdna_kallisto_index/")
abundance_gene_file=gsub( "_transcripts", "_genes",abundance_transcripts_file)

if (!file.exists(tx2gene_file)) {
	stop("tx2gene_file does not exist!")
} else {
	tx2gene=read.delim(tx2gene_file, as.is=T)
	txi.abundance <- tximport(abundance_transcripts_file, type = "kallisto", tx2gene = tx2gene)
	txi.abundance=as.data.frame(txi.abundance)
	txi.abundance=cbind(gene_id=rownames(txi.abundance), txi.abundance)
	write.table(txi.abundance, abundance_gene_file, sep="\t", col.names=T, row.names=F, quote=F)
}


