# Performs differential gene expression with DESeq2
# Written by Edouard Henrion - May 2017
# Usage : Rscript deseq2.R -d path_design -c path_rawcountfile -o output_dir

library(DESeq2)
library(methods)

# Usage

usage = function(errM) {
    cat("\nUsage : Rscript deseq2.R [option] <Value>\n")
    cat("       -d      : design file\n")
    cat("       -c      : raw count file\n")
    cat("       -b      : batch file\n")
    cat("       -o      : output directory\n")
    cat("       -h      : this help\n\n")

    stop(errM)
}

set.seed(123456789)
perform_dge = function(counts, groups, batch, count_limit, path) {

    # Retain row which have > count_limit
    counts <- round(counts[rowSums(counts) > count_limit, ])

    # Normalize and do test
    coldata = data.frame(row.names=colnames(counts), condition=groups)
    ddsFullCountTable = DESeq2::DESeqDataSetFromMatrix(countData = counts, colData=coldata, design=~condition+batch)

    dds <- DESeq2::DESeq(ddsFullCountTable)
    res <- DESeq2::results(dds)
    res[, 1] = row.names(res)
    res[, 5] = as.numeric(format(res[, 5], digits=2))
    res[, 6] = as.numeric(format(res[, 6], digits=2))
    colnames(res)[c(1, 5, 6)] = c("id", "deseq2.p-value", "deseq2.adj.pvalue")
    write.table(res[order(res[, 6]), c(1, 5, 6)], paste(path, "deseq2_results.csv", sep="/"), quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=TRUE)
    fileOpen = paste(path, "edger_results.csv", sep="/")
    d1 <- read.table(fileOpen, header=T, sep="\t", quote="", comment.char="")
    res <- as.data.frame(res)
    d2 <- merge(d1, res[, c(1, 5, 6)], by.x=1, by.y=1, sep="\t")
    d2 <- d2[order(d2[, (ncol(d2)-1)]), ]
    vecWrite <- c(1:4, (ncol(d2)-1), ncol(d2), 5:6, 7:(ncol(d2)-2))
    write.table(d2[, vecWrite], paste(path, "dge_results.csv", sep="/"), quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=TRUE)
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
batches = ""
if (batch_file != "") {
    batches = read.csv2(batch_file, header=T, sep="\t", na.strings="0", check.names=F, colCalsses=c('character', 'character'))
    # make sure design and batch are following the same sample order
    merge_sorted <- merge(design, batches, sort=F, all.x=T, by.x=1)
    batches <- merge_sorted[, 1:ncol(merge_sorted)]

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

    batch = ""
    if (!(batches == "")) {
        batch = as.character(batches)[!(is.na(current_design))]
    }

    cat("Processing for the design\n")
    cat(paste("Name folder: ", name_folder, "\n", sep=""))
    cat(paste("Design : ", paste(subsampleN, group, sep="=", collapse=" ; "), "\n", spe=""))

    # Perform gene differential expressoin
    perform_dge(current_countMatrix, group, batch, count_limit, name_folder)
}
