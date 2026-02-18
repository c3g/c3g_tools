# Performs differential gene expression with DESeq2
# Written by Edouard Henrion - May 2017
# Usage : Rscript deseq2.R -d path_design -c path_rawcountfile -o output_dir
#
# Modified by Senthil - Feb 2026
# Changes to merge block (lines 53-86):
#   - Prefixed edgeR columns with "edgeR_" (log_FC, log_CPM, p.value, adj.p.value)
#   - Included all DESeq2 result columns (baseMean, log2FoldChange, lfcSE, stat,
#     pvalue, padj) with "deseq2_" prefix in dge_results.csv
#   - Column order: id | gene_symbol | edgeR stats | deseq2 stats | sample counts
#   - Sorted by deseq2_padj ascending
#   - deseq2_results.csv output remains unchanged

library(DESeq2)
library(methods)

# Usage

usage = function(errM) {
    cat("\nUsage : Rscript deseq2.R [option] <Value>\n")
    cat("       -d      : design file\n")
    cat("       -c      : raw count file\n")
    cat("       -b      : batch file\n")
    cat("       -o      : output directory\n")
    cat("       -l      : perform local fit instead of default parametric dispersion fit\n")
    cat("       -h      : this help\n\n")

    stop(errM)
}

set.seed(123456789)
perform_dge = function(counts, groups, batch, count_limit, path, locfit = 0) {

    # Retain row which have > count_limit
    counts <- round(counts[rowSums(counts) > count_limit, ])

    # Normalize and do test
    if ((length(batch)==1) && (batch=="")) {
        coldata = data.frame(row.names=colnames(counts), condition=groups)
        ddsFullCountTable = DESeq2::DESeqDataSetFromMatrix(countData = counts, colData=coldata, design=~condition)
    } else {
        coldata = data.frame(row.names=colnames(counts), condition=groups, batch=batch)
        ddsFullCountTable = DESeq2::DESeqDataSetFromMatrix(countData = counts, colData=coldata, design=~batch+condition)
    }

    if (locfit == 1) {
        cat("Using fitType='local' for dispersion estimation\n")
        dds <- DESeq2::DESeq(ddsFullCountTable, fitType = "local")
    } else {
        dds <- DESeq2::DESeq(ddsFullCountTable)
    }
    res <- DESeq2::results(dds)

    # Convert to data.frame and add id column (do NOT overwrite baseMean)
    res_df <- as.data.frame(res)
    res_df <- cbind(id = row.names(res_df), res_df)
    rownames(res_df) <- NULL

    # Write full DESeq2 results (7 columns: id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
    write.table(res_df[order(res_df$padj), ], paste(path, "deseq2_results.csv", sep="/"), quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=TRUE)

    # Merge with edgeR results for combined dge_results.csv
    if ((length(batch)==1) && (batch=="")) {
        fileOpen = paste(path, "edger_results.csv", sep="/")
        if (file.exists(fileOpen)) {
            d1 <- read.table(fileOpen, header=T, sep="\t", quote="", comment.char="")

            # Prefix edgeR statistical columns (columns 3-6), leave id, gene_symbol, and sample counts as-is
            edger_stat_cols <- c("log_FC", "log_CPM", "edger.p.value", "edger.adj.p.value")
            colnames(d1)[3:6] <- paste0("edgeR_", edger_stat_cols)

            # Prepare full DESeq2 data with prefixed columns
            deseq2_for_merge <- data.frame(
                id = res_df$id,
                deseq2_baseMean = res_df$baseMean,
                deseq2_log2FoldChange = res_df$log2FoldChange,
                deseq2_lfcSE = res_df$lfcSE,
                deseq2_stat = res_df$stat,
                deseq2_pvalue = res_df$pvalue,
                deseq2_padj = res_df$padj
            )

            d2 <- merge(d1, deseq2_for_merge, by.x=1, by.y=1)

            # Column order: id | gene_symbol | edgeR stats | deseq2 stats | sample counts
            n_edger <- ncol(d1)
            edger_stat_idx <- 2:6                        # gene_symbol + 4 edgeR stat cols
            deseq2_stat_idx <- (n_edger + 1):(n_edger + 6)  # 6 deseq2 columns appended by merge
            sample_idx <- 7:n_edger                      # sample count columns from edgeR table
            d2 <- d2[order(d2$deseq2_padj), ]
            write.table(d2[, c(1, edger_stat_idx, deseq2_stat_idx, sample_idx)], paste(path, "dge_results.csv", sep="/"), quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=TRUE)
        } else {
            cat(paste("Warning: edgeR results not found at", fileOpen, "- skipping dge_results.csv merge\n"))
        }
    }
}

##################################

ARG = commandArgs(trailingOnly = T)

## default arg values
count_limit = 9
design_file = ""
batch_file = ""
rawcount_file = ""
out_path = ""
locfit = 0

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-d") {
        design_file = ARG[i+1]
    } else if (ARG[i] == "-b") {
        batch_file = ARG[i+1]
    } else if (ARG[i] == "-c") {
        rawcount_file = ARG[i+1]
    } else if (ARG[i] == "-o") {
        out_path = ARG[i+1]
    } else if (ARG[i] == "-l") {
        locfit = 1
    } else if (ARG[i] == "-h") {
        usage("")
    }
}

## check arg consitency
if (!(file.exists(design_file))) {
    usage("Error : Design file not found")
}
if (!(file.exists(rawcount_file))) {
    usage("Error : Raw count file not found")
}
if (out_path == "") {
    usage("Error : Output directory not specified")
}
if (batch_file != "") {
    if (!(file.exists(batch_file))) {
        usage("Error : Batch file not found")
    }
}

# remove trailing "/" if necessary
tmpOP = strsplit(out_path, "")
if (tmpOP[[1]][length(tmpOP[[1]])] == "/") {
    out_path = paste(tmpOP[[1]][1:(length(tmpOP[[1]])-1)], collapse="")
}

design = read.csv2(design_file, header=T, sep="\t", na.strings="0", check.names=F, colClasses=c('character', rep('numeric', unique(count.fields(design_file))-1)))
rawcount = read.csv(rawcount_file, header=T, sep="\t", check.names=F)

print(design)

name_sample = as.character(as.vector(design[, 1]))
countMatrix = round(rawcount[, 3:ncol(rawcount)])

# Check if a batch file has to be used
batch = ""
if (batch_file != "") {
    batches = read.csv2(batch_file, header=T, sep="\t", na.strings="0", check.names=F, colClasses=c('character', 'character'))
    # make sure design and batch are following the same sample order
    merge_sorted = merge(design, batches, sort=F, all.x=T, by.x=1)
    batch = merge_sorted[ncol(merge_sorted)]
    print(batch)
}

# Iterate over each design
for (i in 2:ncol(design)) {

    name_folder = paste(out_path, names(design[i]), sep="/")

    # Create output directory
    if (!file.exists(name_folder)) {
        system(paste("mkdir", name_folder, sep=" "))
    }

    current_design = design[, i]
    subsampleN = name_sample[!(is.na(current_design))]
    group = as.character(current_design)[!(is.na(current_design))]
    current_countMatrix = NULL
    for (j in 1:length(subsampleN)) {
        current_countMatrix = cbind(current_countMatrix, countMatrix[, is.element(colnames(countMatrix), subsampleN[j])])
    }
    colnames(current_countMatrix) = subsampleN
    rownames(current_countMatrix) = rawcount[, 1]

    cat("Processing for the design\n")
    cat(paste("Name folder: ", name_folder, "\n", sep=""))
    cat(paste("Design : ", paste(subsampleN, group, sep="=", collapse=" ; "), "\n", spe=""))

    # Perform gene differential expressoin
    perform_dge(current_countMatrix, group, batch, count_limit, name_folder, locfit)
}