# Performs differential transcript expression with Ballgown 
# Compatible only with stringtie results. 
# Written by Hector Galvez
# Usage: Rscript ballgown.R -d path_design -o output_dir -g gtf
# R module versions: mugqic/R_Bioconductor/3.5.0_3.7 or mugqic/R_Bioconductor/3.5.1_3.7

library(ballgown)
library(methods)
library(dplyr)
library(tibble)

# Usage 

usage=function(errM) { 
    cat("\nUsage : Rscript ballgown.R [option] <Value>\n")
    cat("       -d       : design file \n")
    cat("       -o       : output directory \n")
    cat("       -h       : this help message\n\n") 
    stop(errM)
} 

set.seed(123456789) 

contrast_oper=function(d, current_design) {
    # SAVE EXPRESSION DATA (IN AVERAGE PER-BASE COVERAGE [cov] AND FPKM) 
    smpl_num = length(d$contrast)
    # Define output names
    trx_cov = paste("ballgown",current_design, paste("transcript", "cov_table", sep="."), sep="/")
    trx_fpkm = paste("ballgown",current_design, paste("transcript", "fpkm_table", sep="."), sep="/")
    trx_all = paste("ballgown",current_design, paste("transcript", "attr_table", sep="."), sep="/")
    gene_all = paste("ballgown",current_design, paste("gene", "attr_table", sep="."), sep="/")
    # Prepare ballgown object
    bg <- ballgown(d$path, meas="all")
    bg.idx <- indexes(bg)
    pData(bg) <- data.frame(id=d$sample, group=d$contrast)
    
    # Generate an index for trx ID, gene ID, and gene name
    gnNames <- data.frame(geneNames(bg))
    gnNames <- add_column(gnNames, t_id = row.names(gnNames), .before = 1) 
    gnNames <- add_column(gnNames, g_id = geneIDs(bg), .after = 1) 
    colnames(gnNames) <- c("t_id", "gene_id", "gene_name")
    
    # Produce expression data tables
    trx_tbl = texpr(bg, 'all') 
    trx_tbl$t_id <- as.character(trx_tbl$t_id)
    gene_tbl = data.frame(gexpr(bg)) 
    gene_tbl = add_column(gene_tbl, gene_id = row.names(gene_tbl), .before = 1)     

    # Save transcript expression data
    write.table(dplyr::select(trx_tbl, t_id:gene_name, contains('cov.')), file=trx_cov, row.names=FALSE, quote=FALSE, sep="\t") 
    write.table(dplyr::select(trx_tbl, t_id:gene_name, contains('FPKM.')), file=trx_fpkm, row.names=FALSE, quote=FALSE, sep="\t") 
    write.table(texpr(bg, 'all'), file=trx_all, row.names=FALSE, quote=FALSE, sep="\t") 
    # Save gene expression data
    write.table(gene_tbl, file=gene_all, row.names=FALSE, quote=FALSE, sep="\t") 
    
    # PERFORM DIFFERENTIAL ANALYSIS
    # Define output names
    dge = paste("ballgown", current_design, "gene_exp.diff", sep="/") 
    dte = paste("ballgown", current_design, "transcript_exp.diff", sep="/") 

    # Perform differential analysis
    dte_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group', getFC = TRUE) 
    colnames(dte_results)[2] <- "t_id"
    dte_results$t_id <- as.character(dte_results$t_id)
    out_dte <- dplyr::left_join(dte_results, dplyr::select(trx_tbl, t_id, gene_id, gene_name, contains('FPKM.')), by="t_id")

    dge_results = stattest(bg, feature='gene', meas='FPKM', covariate='group', getFC = TRUE) 
    colnames(dge_results)[2] <- "gene_id"
    dge_results$gene_id <- as.character(dge_results$gene_id)
    out_dge <- dplyr::left_join(dge_results, gene_tbl, by="gene_id")

    # Write outputs
    write.table(dplyr::arrange(dplyr::select(out_dte, -feature), qval), file=dte, row.names=FALSE, quote=FALSE, sep="\t")
    write.table(dplyr::arrange(dplyr::select(out_dge, -feature), qval), file=dge, row.names=FALSE, quote=FALSE, sep="\t")
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
}s

# Read in the design file
design = read.csv2(design_file, header=T, sep="\t", na.strings= "0", check.names=F, colClasses = c('character', rep('numeric',unique(count.fields(design_file))-1)))

name_sample = as.character(as.vector(design[,1])) 

for (i in 2:ncol(design)) {

    name_folder <- paste(out_path, names(design[i]), sep="/")

    if (!file.exists(name_folder)) {
        dir.create(name_folder, showWarnings=F, recursive=T)
    }

    current_design <- design[,i]
    subsampleN <- name_sample[!(is.na(current_design))]
    group <- as.character(current_design)[!(is.na(current_design))]
    stringtiePaths <- paste("stringtie", subsampleN, sep = "/")
    groupN <- unique(group)
    
    # Define current sample to contrast (s2c) table, which will be the main input for dte and dge functions
    current_s2c <- tibble(sample = subsampleN, contrast = group , path = stringtiePaths)
    print(current_s2c)

    cat("Processing for the design\n")
    cat(paste("Name folder: ", name_folder, "\n", sep=""))
    cat(paste("Design: ", paste(subsampleN, group,sep="=", collapse=" ; "), "\n", sep=""))

    contrast_oper(current_s2c, names(design[i]))
}


