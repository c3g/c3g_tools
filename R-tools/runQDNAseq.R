# runs QDNAseq analysis for chromosomal aberrations
# Author: Robert Eveleigh, 2025
# Usage: Rscript runQDNAseq.R -i <input_bam> -o <out_dir> -b <int> -r <hg38> -s <sample_name>

# Load necessary libraries
library(QDNAseq)

# Usage
usage = function(errM) {
  cat("\nUsage : Rscript -i <input> -o <out_dir> -s <sample_name> [options]\n")
  cat("     -i  : bam file\n")
  cat("     -o  : output directory\n")
  cat("     -b  : bin size, default = 15\n")
  cat("     -r  : reference to use, default = hg19\n")
  cat("     -s  : sample name\n\n")

  stop(errM)
}

ARG <- commandArgs(trailingOnly = T)

# default arg values
bamfile <- ""
outdir <- ""
binsize <- 15
reference <- "hg19"
sample <- ""

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-i") {
        bamfile <- ARG[i+1]
    } else if (ARG[i] == "-o") {
        outdir <- ARG[i+1]
    } else if (ARG[i] == "-b") {
        binsize <- ARG[i+1]
    } else if (ARG[i] == "-r") {
        reference <- ARG[i+1]
    } else if (ARG[i] == "-s") {
        sample <- ARG[i+1]
    } else if (ARG[i] == "-h") {
        usage("")
    }
}

## check arg consitency
if (!(file.exists(bamfile))) {
    usage("Error : Bam file not found")
}
if (out_path == "") {
    usage("Error : Output directory not specified")
}
if (sample == "") {
    usage("Error : Sample name not specified")
}

# Check if the output directory exists, if not create it.
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

library(paste("QDNAseq.", args$reference, sep = ""))

# STEP 1: Load bins with specified size
if (reference == "hg19") {
  bins <- getBinAnnotations(binSize = binsize)
} else if (reference == "hg38") {
  bins <- getBinAnnotations(binSize = binsize, genome = "hg38")
} else {
  bins <- getBinAnnotations(binSize = binsize, genome = reference)
}

# STEP 2: Process BAM files
options(future.globals.maxSize = 2000 * 1024^2) # Set to 2GB (adjust as needed)

readCounts <- binReadCounts(bins, bamfiles = bamfile, chunkSize = "chr")

# Helper function to generate a *single* output file name
make_output_filename <- function(basename, format = "pdf") {
    file.path(outdir, paste0(sample, ".", basename, ".", binsize, "k.", reference, ".", format))
}

pdf(make_output_filename("read_counts_per_bin"), width =  8, height = 4)
plot(readCounts, logTransform = FALSE, ylim = c(-50, 2000))
dev.off()


readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)

pdf(make_output_filename("isobarPlot"), width = 8, height = 4)
isobarPlot(readCountsFiltered)
dev.off()

readCountsFiltered <- estimateCorrection(readCountsFiltered)

pdf(make_output_filename("noise_plot"), width = 8, height = 4)
noisePlot(readCountsFiltered)
dev.off()

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

pdf(make_output_filename("copy_number_smooth"), width = 8, height = 4)
plot(copyNumbersSmooth)
dev.off()

# Exporting - loop through formats if needed
formats_to_export <- c("tsv", "igv", "bed")
for (format in formats_to_export) {
  exportBins(copyNumbersSmooth, file = make_output_filename("CNV", format), format = format)
}

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

pdf(make_output_filename("copy_number_segmented"), width = 8, height = 4)
plot(copyNumbersSegmented)
dev.off()

copyNumbersCalled <- callBins(copyNumbersSegmented)

pdf(make_output_filename("copy_number_calls"), width = 8, height = 4)
plot(copyNumbersCalled)
dev.off()

# Exporting - loop through formats if needed
formats_to_export <- c("vcf", "seg")
for (format in formats_to_export) {
  exportBins(copyNumbersCalled, file = make_output_filename("CNV_calls", format), format = format)
}
