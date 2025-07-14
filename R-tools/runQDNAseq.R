#!/usr/bin/env Rscript

# Load necessary libraries
library(QDNAseq)
library(QDNAseq.hg38)
library(argparse)

# Set up argument parsing
parser <- ArgumentParser(description = "Perform CNV analysis using QDNAseq")

parser$add_argument("-i", "--input_bam", dest="bamfile", required=TRUE, help="Path to the input BAM file")
parser$add_argument("-o", "--outdir", dest="outdir", required=TRUE, help="Path to the output directory")
parser$add_argument("-b", "--binsize", dest="binsize", type="integer", default=15, help="Bin size in kb (default: 15)") # Added binsize argument
parser$add_argument("-r", "--reference", dest="reference", default="hg19", choices=c("hg19", "hg38"), help="Reference genome (hg19 or hg38, default: hg19)") # Added reference genome argument
parser$add_argument("-s", "--sample", dest="sample", required = TRUE, help="Name of the sample, used to name output files")

args <- parser$parse_args()

# Check if the output directory exists, if not create it.
if (!dir.exists(args$outdir)) {
    dir.create(args$outdir, recursive = TRUE) # recursive = TRUE will create parent folders if needed
}

# STEP 1: Load bins with specified size
if (args$reference == "hg19") {
  bins <- getBinAnnotations(binSize = args$binsize)
} else if (args$reference == "hg38") {
  bins <- getBinAnnotations(binSize = args$binsize, genome = "hg38") # Specify reference genome
}

# STEP 2: Process BAM files 
options(future.globals.maxSize = 2000 * 1024^2) # Set to 2GB (adjust as needed)

readCounts <- binReadCounts(bins, bamfiles = args$bamfile, chunkSize = "chr")

# Helper function to generate a *single* output file name
make_output_filename <- function(basename, format = "pdf") {
    file.path(args$outdir, paste0(sample, ".", basename, ".", args$binsize, "k.", args$reference, ".", format))
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
