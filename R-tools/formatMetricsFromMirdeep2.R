# Format and report properly as a table all the metrics of a mirdeep2 analysis
# Written by Edouard Henrion - March 2017
# Usage : Rscript formatMetricsFromMirdeep2.R -raw raw_read_cmd -trim post_trim_cmd -mapg map_genome_cmd -filter filter_reads_cmd -mapm map_mirna_cmd -format r_format_cmd

options(stringsAsFactors=F)
library(gdata)
library(magrittr)

usage=function(errM) {
    cat("\nUsage : Rscript mirdeep2Metrics.R [option] <Value>\n")
    cat("       -raw      : list of arf_files\n")
    cat("       -trim     : list of fasta_files\n")
    cat("       -mapg     :
    cat("       -filter   :
    cat("       -mapm     :
    cat("       -format   :
    cat("       -h        : this help\n\n")

    stop(errM)
}

ARG = commandArgs(trailingOnly = T)
## default arg values
target_species=""
other_species=""
mirbase_path=""
## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-raw") {
        raw_read_cmd=ARG[i+1]
    } else if (ARG[i] == "-trim") {
        post_trim_cmd=ARG[i+1]
    } else if (ARG[i] == "-mapg") {
        map_genome_cmd=ARG[i+1]
    } else if (ARG[i] == "-filter") {
        filter_reads_cmd=ARG[i+1]
    } else if (ARG[i] == "-mapm") {
        map_mirna_cmd=ARG[i+1]
    } else if (ARG[i] == "-format") {
        r_format_cmd=ARG[i+1]
    } else if (ARG[i] == "-h") {
        usage("")
    }
}
## check arg consitency
if (raw_read_cmd == "") {
    usage("Error : Output directory not specified")
}
if (post_trim_cmd == "") {
    usage("Error : Output directory not specified")
}
if (map_genome_cmd == "") {
    usage("Error : Output directory not specified")
}
if (filter_reads_cmd == "") {
    usage("Error : Output directory not specified")
}
if (map_mirna_cmd == "") {
    usage("Error : Output directory not specified")
}
if (r_format_cmd == "") {
    usage("Error : Output directory not specified")
}


# This function parses output of wc -l
parse.linecount <- function(filename, div.by=1)
{
x = readLines(filename) %>% gdata::trim(.) %>% strsplit(split=' ') %>% do.call(rbind,.) %>% as.data.frame
x = x[-nrow(x),]
r = as.numeric(x[,1]); names(r) = x[,2]
r / div.by
}

# this function sums up over mirdeep2 style headers and translates codes
parse.mirdeep2.header <- function(filename)
{
x = readLines(filename)
x = unique(x) # multimapped
x = data.frame("code" = gsub("_.*", "", x), "count" = as.numeric(gsub(".*_x", "", x)) )
x = aggregate(x[,2], list(x[,1]), sum)
r = x[,2]; names(r) = x[,1]
r
}

codes = read.delim("readsets_codes.tsv", header=F)
rownames(codes) = codes[,2]

raw = parse.linecount("metrics/raw_reads_fastq.linecount", 4)
names(raw) = basename(dirname(names(raw)))

trimmed =  parse.linecount("metrics/reads_trimmed_fastq.linecount", 4)
names(trimmed) = basename(dirname(names(trimmed)))

length.filtered = parse.mirdeep2.header("metrics/length_filtered_fasta.header")
names(length.filtered) = codes[names(length.filtered),1]

mapped = parse.mirdeep2.header("metrics/mapped_arf.header")
names(mapped) = codes[names(mapped),1]

mirna.mapped = parse.mirdeep2.header("metrics/miRNA_aligned_arf.header")
names(mirna.mapped) = codes[names(mirna.mapped),1]

counts = do.call(cbind, list(raw[codes[,1]], trimmed[codes[,1]], length.filtered[codes[,1]], mirna.mapped[codes[,1]], mapped[codes[,1]]))
colnames(counts) = c("Raw", "Trimmed and Adapter Filtered", "Minimum Length Filtered", "miRNA Mapped", "Genome Mapped")
save(counts, file="metrics/read_counts.RData")
write.csv(format(counts, big.mark=','), file="metrics/read_counts.csv")
