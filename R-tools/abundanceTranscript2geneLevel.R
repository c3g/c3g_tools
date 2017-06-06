#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

abundance_transcripts_file=args[1]
gtf_file=args[2]


# abundance_transcripts_file="/Users/emercier/abacus/projects/RNAseq_light_BFXDEV_615/RNAseq_light2/Pipeline_output/kallisto/KW390ASoy200/abundance_transcripts.tsv"
# gtf_file="/cvmfs/soft.mugqic/CentOS6/genomes/species/Mus_musculus.GRCm38/annotations/Mus_musculus.GRCm38.Ensembl83.transcript_id.gtf"
require(tximport)

index_path=file.path(dirname(gtf_file),"cdna_kallisto_index/")
tx2gene_file_repo=paste0(index_path, basename(tools::file_path_sans_ext(gtf_file)), ".tx2gene")
tx2gene_file_local=paste0(dirname(dirname(abundance_transcripts_file)),"/", basename(tools::file_path_sans_ext(gtf_file)), ".tx2gene") # ../kallisto
abundance_gene_file=gsub( "_transcripts", "_genes",abundance_transcripts_file)

if (file.exists(tx2gene_file_repo)) {
	tx2gene_file=tx2gene_file_repo
} else if (file.exists(tx2gene_file_local)) {
	tx2gene_file=tx2gene_file_local
} else {
	require(rtracklayer)
	print("Building transcripts2genes...")

	gtf=import(gtf_file, format = "gff2")

	tx2gene=cbind(tx_id=gtf$transcript_id, gene_id=gtf$gene_id) #gene_name
	tx2gene=tx2gene[!is.na(tx2gene[,1]),]
	tx2gene=unique(tx2gene)
	tx2gene=as.data.frame(tx2gene)

	write.table(x=tx2gene, file=tx2gene_file_local, sep="\t", col.names=T, row.names=F, quote=F)
	tx2gene_file=tx2gene_file_local
}


tx2gene=read.delim(tx2gene_file, as.is=T)
txi.abundance <- tximport(abundance_transcripts_file, type = "kallisto", tx2gene = tx2gene)
txi.abundance=as.data.frame(txi.abundance)
txi.abundance=cbind(gene_id=rownames(txi.abundance), txi.abundance)
write.table(txi.abundance, abundance_gene_file, sep="\t", col.names=F, row.names=T, quote=F)
