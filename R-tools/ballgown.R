# Performs differential transcript expression with Ballgown 
# Compatible only with stringtie results. 
# Written by Hector Galvez
# Usage: Rscript ballgown.R -d path_design -o output_dir -g gtf
# R module versions: mugqic/R_Bioconductor/3.5.0_3.7 or mugqic/R_Bioconductor/3.5.1_3.7

library(ballgown)
library(methods)
library(dplyr)
library(tibble)
library(magrittr)
library(readr)

# Usage 
usage=function(errM) { 
    cat("\nUsage : Rscript ballgown.R [option] <Value>\n")
    cat("       -d       : design file \n")
    cat("       -o       : output directory \n")
    cat("       -h       : this help message\n\n") 
    stop(errM)
} 

set.seed(123456789) 

contrast_oper=function(d, current_design, out_folder) {
    # SAVE EXPRESSION DATA (IN AVERAGE PER-BASE COVERAGE [cov] AND FPKM) 
    smpl_num = length(d$contrast)
    ## Define output names
    trx_cov = file.path(out_folder, paste("transcript", "cov_table", sep="."))
    trx_fpkm = file.path(out_folder, paste("transcript", "fpkm_table", sep="."))
    trx_all = file.path(out_folder, paste("transcript", "attr_table", sep="."))
    gene_all = file.path(out_folder, paste("gene", "attr_table", sep="."))

    ## TODO: Double-check that the files actually exist at d$path. Give helpful error if files can't be found.
    ## Prepare ballgown object
    bg <- ballgown(d$path, meas="all")
    bg.idx <- indexes(bg)
    pData(bg) <- data.frame(id=d$sample, group=d$contrast)
    
    ## Generate an index for trx ID, gene ID, and gene name
    gnNames <- data.frame(geneNames(bg))
    gnNames <- add_column(gnNames, t_id = row.names(gnNames), .before = 1) 
    gnNames <- add_column(gnNames, g_id = geneIDs(bg), .after = 1) 
    colnames(gnNames) <- c("t_id", "gene_id", "gene_name")
    
    ## Produce expression data tables
    trx_tbl = texpr(bg, 'all') 
    trx_tbl$t_id <- as.character(trx_tbl$t_id)
    gene_tbl = data.frame(gexpr(bg)) 
    gene_tbl = add_column(gene_tbl, gene_id = row.names(gene_tbl), .before = 1)     

    ## Save transcript expression data
    trx_tbl %>%
        dplyr::select(t_id:gene_name, contains('cov.')) %>%
        write_tsv(trx_cov)
    
    trx_tbl %>%
        dplyr::select(t_id:gene_name, contains('FPKM.')) %>%
        write_tsv(trx_fpkm)
    
     write_tsv(texpr(bg, 'all'), trx_all)

    ## Save gene expression data
    write_tsv(gene_tbl, path = gene_all)

    # PERFORM DIFFERENTIAL ANALYSIS
    ## Perform differential analysis
    stattest(bg, feature='transcript', meas='FPKM', covariate='group', getFC = TRUE) %>%
        dplyr::mutate(t_id = as.character(id)) %>%
        select(-id) -> dte_results
    
    stattest(bg, feature='gene', meas='FPKM', covariate='group', getFC = TRUE) %>%
        dplyr::mutate(gene_id = as.character(id)) %>%
        select(-id) %>%
        dplyr::left_join(gene_tbl, by="gene_id") -> out_dge
    
    ## Define output names
    dge = file.path(out_folder, "gene_exp.diff") 
    dte = file.path(out_folder, "transcript_exp.diff") 
    
    ## Write outputs
    trx_tbl %>%
        dplyr::select(t_id, gene_id, gene_name, contains('FPKM.')) %>%
        dplyr::left_join(dte_results, ., by="t_id") %>%
        dplyr::select(-feature) %>%
        dplyr::arrange(qval) %>%
        write_tsv(dte)
    
    out_dge %>%
        dplyr::select(-feature) %>%
        dplyr::arrange(qval) %>%
        write_tsv(dge)
}

#########################

ARG = commandArgs(trailingOnly = T) 
## default arg values 
fpath="."
design_file=""
out_path=""
## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-d") {
        design_file=ARG[i+1]
    } else if (ARG[i] == "-o") {
        out_path=ARG[i+1]
    } else if (ARG[i] == "-h") {
        usage("")
    }
}

# Read in the design file
design = read.csv2(design_file, header=T, sep="\t", na.strings= "0", check.names=F, colClasses = c('character', rep('numeric',unique(count.fields(design_file, sep = "\t"))-1)))

name_sample = as.character(as.vector(design[,1])) 

for (i in 2:ncol(design)) {

    name_folder <- file.path(out_path, names(design[i]))

    if (!file.exists(name_folder)) {
        dir.create(name_folder, showWarnings=F, recursive=T)
    }

    current_design <- design[,i]
    subsampleN <- name_sample[!(is.na(current_design))]
    group <- as.character(current_design)[!(is.na(current_design))]
    stringtiePaths <- file.path("stringtie", subsampleN)
    groupN <- unique(group)
    
    # Define current sample to contrast (s2c) table, which will be the main input for dte and dge functions
    current_s2c <- tibble(sample = subsampleN, contrast = group , path = stringtiePaths)
    print(current_s2c)

    cat("Processing for the design\n")
    cat(paste("Output folder: ", name_folder, "\n", sep=""))
    cat(paste("Design: ", paste(subsampleN, group,sep="=", collapse=" ; "), "\n", sep=""))

    contrast_oper(current_s2c, names(design[i]), name_folder)
}


