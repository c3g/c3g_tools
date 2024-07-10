#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

report_dir=args[1]
output_file=args[2]

print("Creating directories and copying outputs...\n")
# Directories and cp data
unlink(file.path(report_dir,"kallisto") , recursive = T)
dir.create(file.path(report_dir,"kallisto"), showWarnings=F,recursive=T)
file.copy(from = "kallisto", to = report_dir, overwrite = T, recursive = T)

library(jsonlite)
options(stringsAsFactors=F)

# Create summary table
print("Creating kallisto summary table for all readsets...\n")
all_readsets_abundance.transcripts <- read.delim(file.path(report_dir, "kallisto/All_readsets/all_readsets.abundance_transcripts.csv"), 
                                                 as.is=T, 
                                                 check.names = FALSE) #to allow the use of dash in sample names

all_readsets_abundance.genes <- read.delim(file.path(report_dir,"kallisto/All_readsets/all_readsets.abundance_genes.csv"), 
                                           as.is=T, 
                                           check.names = FALSE)

all_readsets_trimming_metrics <- read.delim("metrics/trimSampleTable.tsv", as.is=T)
colnames(all_readsets_trimming_metrics) <- c("Sample","Raw Reads #", "Surviving Reads #","Surviving Reads %")
rownames(all_readsets_trimming_metrics) <- all_readsets_trimming_metrics$Sample
readset_names <- all_readsets_trimming_metrics$Sample

mat_colnames <- c("Sample","RawReads", "SurvivingReads", "PercSurviving", "Transcriptome", "Transcripts", "TranscriptsReads","PercTranscriptsReads", "Genes", "GenesReads","PercGenesReads")
mat <- data.frame(matrix(nrow=length(readset_names), ncol=length(mat_colnames)))
rownames(mat) <- readset_names
colnames(mat) <- mat_colnames
mat[,1] <- readset_names

for (i in seq(length(readset_names))){
    readset <- readset_names[i]

    run_info <- fromJSON(file.path(report_dir, "kallisto", readset, "run_info.json"))
    if (grepl("--single", run_info$call)) {
        n_fastq = 1
    } else {
        n_fastq = 2
    } #check if paired or single
    
    mat[i,"RawReads"] <- all_readsets_trimming_metrics[readset, "Raw Reads #"]
    mat[i,"SurvivingReads"] <- all_readsets_trimming_metrics[readset,"Surviving Reads #"]
    mat[i,"PercSurviving"] <- signif(all_readsets_trimming_metrics[readset,"Surviving Reads %"],3)

    mat[i,"Transcriptome"] <- run_info$n_targets

    transcript_count <- all_readsets_abundance.transcripts[,readset]
    mat[i,"Transcripts"] <- sum(transcript_count>=5)
    mat[i, "TranscriptsReads"] <- sum(transcript_count)*n_fastq #multiply by 2 if paired
    mat[i,"PercTranscriptsReads"] <- signif(mat[i,"TranscriptsReads"]/mat[i,"RawReads"],3)*100

    gene_count <- all_readsets_abundance.genes[,readset]
    mat[i,"Genes"] <- sum(gene_count>=5)
    mat[i,"GenesReads"] <- sum(gene_count)*n_fastq #multiply by 2 if paired
    mat[i, "PercGenesReads"] <- signif(mat[i,"GenesReads"]/mat[i,"RawReads"],3)*100
}

print("Writing summary table to file...")
write.table(mat, file = output_file, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
