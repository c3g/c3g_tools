## running ggplots2 on hap.py output
## Robert Eveleigh - robert.eveleigh@mcgill.ca
## 2021/06/09
## usage:
##   Rscript run_annofuse.R <ARRIBA_INPUT> <STAR_FUSION_INPUT> <SAMPLE_NAME>

##get args
args=commandArgs(TRUE)
arriba_input=args[1]
star_fusion_input=args[2]
sample_name=args[3]

outfile = paste(sample_name, "putative_driver_fusions.tsv", sep=".")

##Libraries
library("annoFuse")
suppressPackageStartupMessages(library("readr"))

fusion_calls<-annoFuse_single_sample(
    arriba_input,
    star_fusion_input,
    expressionFile = NULL,
    expressionFilter = 1,
    tumorID = sample_name,
    readingFrameFilter = "in-frame|frameshift|other",
    readthroughFilter = FALSE,
    artifactFilter = "GTEx_recurrent_StarF2019|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
    junctionReadCountFilter = 1,
    spanningFragCountFilter = 10
)

# Add reference gene list containing known oncogenes, tumor suppressors, kinases, and transcription factors
geneListReferenceDataTab <- read.delim(
  system.file("extdata", "genelistreference.txt", package = "annoFuseData"),
  stringsAsFactors = FALSE
)
# Add fusion list containing previously reported oncogenic fusions.
fusionReferenceDataTab <- read.delim(
  system.file("extdata", "fusionreference.txt", package = "annoFuseData"),
  stringsAsFactors = FALSE
)

putative_driver_fusions <- fusion_driver(
  standardFusioncalls = fusion_calls, 
  annotated = TRUE, 
  geneListReferenceDataTab = geneListReferenceDataTab, 
  fusionReferenceDataTab = fusionReferenceDataTab,
  checkDomainStatus = FALSE
)

putative_driver_fusions <- aggregate_fusion_calls(
    putative_driver_fusions, 
    removeother = FALSE
)


write_tsv(putative_driver_fusions, outfile, append=TRUE, col_names=TRUE)
