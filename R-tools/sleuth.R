# Performs differential transcript and gene expression with Sleuth
# Compatible only with kallisto results. 
# Written by Hector Galvez
# Usage: Rscript sleuth.R -d path_design -o output_dir -t tx2gene

library(readr)
library(sleuth)
library(cowplot)
library(methods)
library(pheatmap)
library(RColorBrewer)

# Usage 

usage=function(errM) { 
    cat("\nUsage : Rscript sleuth.R [option] <Value>\n")
    cat("       -d       : design file\n") 
    cat("       -o       : output directory\n")
    cat("       -t       : transcript to gene table\n")
    cat("       -h       : this help message\n\n") 
    stop(errM)
} 

set.seed(123456789) 

perform_dte=function(d, current_design) {
    # Define output names
    lrt_out = paste("sleuth",current_design, paste("results", "lrt", "trx", "csv", sep="."), sep="/")
    wt_out = paste("sleuth",current_design, paste("results", "wt", "trx", "csv", sep="."), sep="/")
    pca_out = paste("sleuth",current_design, paste("pca_plot", "trx", "png", sep="."), sep="/") 
    # Prepare sleuth object. Use a different transformation function so that we can compute log2 fold change. 
    so <- sleuth_prep(d, ~ contrast, transformation_function = function(x) log2(x + 0.5), num_cores=2)
    # Perform likelihood ratio test (LRT)
    so <- sleuth_fit(so, ~ contrast, "full") 
    so <- sleuth_fit(so, ~ 1, "reduced") 
    so <- sleuth_lrt(so, "reduced", "full") 
    lrt_results <- sleuth_results(so, "reduced:full", "lrt", show_all=TRUE) 
    lrt_results <- merge(lrt_results, ttg, by = c("target_id"), all.x = TRUE)
    # Perform wald test (WT)
    so <- sleuth_wt(so, "contrast2", "full")
    wt_results <- sleuth_results(so, "contrast2", "wt", show_all=TRUE) 
    wt_results <- merge(wt_results, ttg, by = c("target_id"), all.x = TRUE) 
        # Add log-change columng to wald test results 
    wt_results$fold_change <- log2(2^wt_results$b)

    # Write results files and plots
    write.csv(lrt_results[order(lrt_results$qval),], file=lrt_out, row.names = FALSE) 
    write.csv(wt_results[order(wt_results$qval),], file=wt_out, row.names = FALSE)
    #ggsave(pca_out, plot_pca(so, color_by = "contrast", text_labels = TRUE), dpi = 300) 
}

perform_dge=function(d, current_design) {
    # Define output names 
    lrt_out = paste("sleuth", current_design, paste("results", "lrt", "gene", "csv", sep="."), sep="/")
    wt_out = paste("sleuth", current_design, paste("results", "wt", "gene", "csv", sep="."), sep="/")
    pca_out = paste("sleuth", current_design, paste("pca_plot", "gene", "png", sep="."), sep="/") 
    htmap_out = paste("sleuth", current_design, paste("heatmap", "topFCgenes", "pdf", sep="."), sep="/") 
    # Prepare sleuth object. 
    so <- sleuth_prep(d, ~ contrast, target_mapping = ttg, aggregation_column = 'ens_gene', transformation_function = function(x) log2(x + 0.5), num_cores=2)
    # Perform likelihood ratio test (LRT) 
    so <- sleuth_fit(so, ~ contrast, "full") 
    so <- sleuth_fit(so, ~ 1, "reduced") 
    so <- sleuth_lrt(so, "reduced", "full")
    lrt_results <- sleuth_results(so, "reduced:full", "lrt", show_all = TRUE) 
    # Perform wald test (WT) 
    so <- sleuth_wt(so, "contrast2", "full") 
    wt_results <- sleuth_results(so, "contrast2", "wt", show_all=TRUE)
    wt_results$fold_change <- log2(2^wt_results$b)

    # Find top and bottom 20 genes by fold-change
    wt_uniq_results <- wt_results[!duplicated(wt_results$target_id),]
    ordered_by_fc <- wt_uniq_results[order(-abs(wt_uniq_results$fold_change)),c("ext_gene","target_id")]
    top40 <- ordered_by_fc[1:40,]
    heatmap_subset <- so$bs_summary$obs_counts[top40$target_id,]
    row.names(heatmap_subset) <- top40$ext_gene

    # Write results files
    write.csv(lrt_results[!duplicated(lrt_results$target_id), !names(lrt_results) %in% c("transcript_version","degrees_free")], file=lrt_out, row.names = FALSE) 
    write.csv(wt_results[!duplicated(wt_results$target_id), !names(wt_results) %in% c("transcript_version","degrees_free")], file=wt_out, row.names = FALSE) 
    # Save plots
    #png(file=pca_out, height=7, width=7, units="in", res=500)
    #plot_pca(so, color_by = "contrast", text_labels = TRUE)
    #dev.off()
    #
    #pdf(file = heatmap_out)
    #pheatmap(heatmap_subset, color =colors, filename= htmap_out, silent)
    #dev.off()
}


#########################

ARG = commandArgs(trailingOnly = T) 
## default arg values 
fpath="."
design_file=""
tx2gene=""
out_path=""
## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-d") {
        design_file=ARG[i+1]
    } else if (ARG[i] == "-o") {
        out_path=ARG[i+1]
    } else if (ARG[i] == "-t") {
        tx2gene=ARG[i+1]
    } else if (ARG[i] == "-h") {
        usage("")
    }
}

# Read in the design file
design = read.csv2(design_file, header=T, sep="\t", na.strings= "0", check.names=F, colClasses = c('character', rep('numeric',unique(count.fields(design_file))-1)))

# Read in the transcript to gene file
ttg <- read_csv(tx2gene) 

# Prepare blue color scheme for heatmaps 
colors <- colorRampPalette( brewer.pal(9, "Blues"))(255)

name_sample = as.character(as.vector(design[,1])) 

for (i in 2:ncol(design)) {

    name_folder <- paste(out_path, names(design[i]), sep="/")

    if (!file.exists(name_folder)) {
        dir.create(name_folder, showWarnings=F, recursive=T)
    }

    current_design <- design[,i]
    subsampleN <- name_sample[!(is.na(current_design))]
    group <- as.character(current_design)[!(is.na(current_design))]
    kallistoPaths <- paste("kallisto", subsampleN, sep = "/")
    groupN <- unique(group)

    current_s2c <- tibble(sample = subsampleN, contrast = group , path = kallistoPaths)
    print(current_s2c)

    cat("Processing for the design\n")
    cat(paste("Name folder: ", name_folder, "\n", sep=""))
    cat(paste("Design: ", paste(subsampleN, group,sep="=", collapse=" ; "), "\n", sep=""))

    perform_dte(current_s2c, names(design[i]))
    perform_dge(current_s2c, names(design[i]))
}


