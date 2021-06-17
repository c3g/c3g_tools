#' ---
#' title: "The report for differential binding analysis comparing `r params$c`"
#' author: ""
#' params:
#'  d: ""
#'  r: ""
#'  o: ""
#'  c: ""
#'  b: "alignment"
#'  p: "peak_call"
#'  minOverlap: 2
#'  dir: "differential_binding"
#'  cur_wd: ""
#'  minMembers: 2
#'  method: "DBA_DESEQ2"
#'  
#' ---



knitr::opts_chunk$set(echo=FALSE)
knitr::opts_chunk$set(message=FALSE)
knitr::opts_knit$set(root.dir = params$cur_wd)


#' # DiffBind Analysis

#first set of codes define YAML parameters. It captures the paramters recived from the outside
# then turn offed echo. Otherwise HTML page has all the codes. If you want to see the code, set these two to true
# then set current working directory to user folder. Otherwise, R script executes where it is stored.

#loading libraries 
library(DiffBind)
library(gdata)
library(DESeq2)
library(data.table)
library(dplyr)

#default help when run using Rscript paramters ...
usage = function(errM) {
  cat("\nUsage : Rscript diffbind.R [option] <Value>\n")
  cat("       -d      : design file\n")
  cat("       -r      : readset file \n")
  cat("       -c      : comparison column name \n")
  cat("       -o      : output file path\n")
  cat("       -b      : bam directory\n")
  cat("       -p      : peak file directory\n")
  cat("       -minMembers      : MinMembers in a group\n")
  cat("       -minOverlap      : minOverlap in a group\n")
  cat("       -h      : this help\n\n")
  
  stop(errM)
}

set.seed(123456789)

#only run when rendering using knit
if(isTRUE(getOption('knitr.in.progress'))){
  design_file = params$d
  readset_file = params$r
  out_path = params$o
  comparison = params$c
  bam_dir = params$b
  peak_dir = params$p
  minoverlap = params$minOverlap
  minmembers = params$minMembers
  out_dir = params$dir
  diff_method = params$method
  
  
  
} else{
#only run when execut using traditional method
##################################

ARG = commandArgs(trailingOnly = T)

## default arg values

design_file = ""
readset_file = ""
out_path = ""
comparison= ""
bam_dir=""
peak_dir=""
minoverlap=2
minmembers=2
diff_method="DBA_DESEQ2"

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-d") {
    design_file = ARG[i+1]
  } else if (ARG[i] == "-r") {
    readset_file = ARG[i+1]
  } else if (ARG[i] == "-o") {
    out_path = ARG[i+1]
  } else if (ARG[i] == "-c") {
    comparison = ARG[i+1]
  } else if (ARG[i] == "-b") {
    bam_dir = ARG[i+1]
  } else if (ARG[i] == "-p") {
    peak_dir = ARG[i+1]
  } else if (ARG[i] == "-h") {
    usage("")
  } else if (ARG[i] == "-minOverlap") {
    minoverlap = ARG[i+1]
  } else if (ARG[i] == "-minMembers") {
    minmembers = ARG[i+1]
  } else if (ARG[i] == "-dir") {
    out_dir = ARG[i+1]
  } else if (ARG[i] == "-method") {
    diff_method = ARG[i+1]
  }
}
}

if (!(file.exists(design_file))) {
  usage("Error : Design file not found")
}
if (!(file.exists(readset_file))) {
  usage("Error : Read set file not found")
}
if (comparison == "") {
  usage("Error : comparison is not specified")
}
if (out_path == "") {
  usage("Error : Output file is not specified")
}
if (!(dir.exists(peak_dir))) {
  usage("Error : directory with macs output files is not found")
}
if (!(dir.exists(bam_dir))) {
  usage("Error : directory with bam files is not found")
}
# remove trailing "/" if necessary
tmpOP = strsplit(out_path, "")
if (tmpOP[[1]][length(tmpOP[[1]])] == "/") {
  out_path = paste(tmpOP[[1]][1:(length(tmpOP[[1]])-1)], collapse="")
}


#Rscript /home/pubudu/projects/rrg-bourqueg-ad/pubudu/chipseq_diff/diff_bind.R -d designfile_chipseq.chr22.txt -r readsets.chipseqTest.kidney.tsv -c H3K4me1_N_T -o differential_binding/diffbind_H3K4me1_N_T_dba.txt

#load design file
design <- fread(design_file)

#load readsets file
readset <- fread(readset_file)

#merge readset and design file
samplesheet <- readset
samplesheet$bamReads <- paste(bam_dir, samplesheet$Sample, samplesheet$MarkName,  paste(samplesheet$Sample, 
                                                                                        samplesheet$MarkName, "sorted.dup.filtered.bam", sep="."), sep="/")

samplesheet$bamControl <- paste(bam_dir, samplesheet$Sample, 'input',  paste(samplesheet$Sample, 
                                                                             "input", "sorted.dup.filtered.bam", sep="."), sep="/")

samplesheet <- select(samplesheet, c("Sample", "MarkName", "bamReads", "bamControl"))
colnames(samplesheet) <- c("Sample", "Factor", "bamReads", "bamControl")

samplesheet$PeakCaller <- "macs"
samplesheet$Peaks <- paste(peak_dir, samplesheet$Sample, samplesheet$Factor,  paste0(samplesheet$Sample, ".", samplesheet$Factor, "_peaks.xls"), sep="/")
samplesheet$PeakFormat <- "macs"


Condition.col <- comparison
#subset by column name passed from a variable
design <- subset(design, design[[Condition.col]]==1 | design[[Condition.col]]==2)
#rename column name passed by a variable
names(design)[names(design) == Condition.col] <- "Condition"


samplesheet <- merge(samplesheet, design, by.x=c("Sample", "Factor"), by.y=c("Sample", "MarkName"))
#' #### Using `r diff_method` as the differential analysis method 
#' #### Sample sheet used for the analysis
print(samplesheet)

#check whether all the files are available
for (i in samplesheet$bamReads) {
  if (!(file.exists(i))) {
    stop(paste("Error : bam file", i, "not found"))
  }
}

for (i in samplesheet$bamControl) {
  if (!(file.exists(i))) {
    stop(paste("Error : control bam file", i, "not found"))
  }
}

for (i in samplesheet$Peaks) {
  if (!(file.exists(i))) {
    
    stop(paste("Error : peak file", i, "not found"))
  } 
}

samplesheet$SampleID <- paste(samplesheet$Sample, samplesheet$Factor, sep="_")
samplesheet$ControlID <- paste(samplesheet$Sample, "input", sep="_")



dba.ob <- dba(sampleSheet=samplesheet, minOverlap=minoverlap)
#' Below table shows information related to samples and macs2 peak files. Such as how many peaks are in each peakset, the total number of unique peaks after merging overlapping ones (in the first line), and the dimensions of the default binding matrix.
dba.ob

#' Using the data from the peak calls, a correlation heatmap can be generated which gives an initial clustering of the samples using the cross-correlations of each row of the binding matrix.
#+ Fig1, fig.cap = "Fig 1: Correlation heatmap, using occupancy (peak caller score) data" , fig.align = "center"
plot(dba.ob)

dba.ob.count <- dba.count(dba.ob, bParallel=F)

#' Next step is to Calculate a binding matrix with scores based on read counts for every sample (affinity scores), rather than confidence scores for only those peaks called in a specific sample (occupancy scores. Two additional columns were added to above table in this step. The first shows the total number of aligned reads for each sample (the "Full" library sizes). The second is labeled FRiP, which stands for Fraction of Reads in Peaks. This is the proportion of reads for that sample that overlap a peak in the consensus peak set, and can be used to indicate which samples show more enrichment overall.
dba.ob.count
info <- dba.show(dba.ob.count)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, PeakReads=round(info$Reads*info$FRiP))

#' For each sample, multiplying the value in theReadscolumn by the correspondingFRiPvaluewill yield the number of reads that overlap a consensus peak
libsizes

#' New correlation heatmap based on the count scores

#+ Fig2, fig.cap = "Fig 2: Correlation heatmap, using affinity (read count) data" , fig.align = "center"
plot(dba.ob.count)

#dba.ob.norm <- dba.normalize(dba.ob.count, method=DBA_DESEQ2)
#print(noquote(diff_method))
if(diff_method=="DBA_DESEQ2"){
dba.ob.norm <- dba.normalize(dba.ob.count, method=DBA_DESEQ2)
} else if (diff_method=="DBA_EDGER") {
  dba.ob.norm <- dba.normalize(dba.ob.count, method=DBA_EDGER)
} else if(diff_method=="DBA_ALL_METHODS"){
  dba.ob.norm <- dba.normalize(dba.ob.count, method=DBA_ALL_METHODS)
}


dba.ob.cont <- dba.contrast(dba.ob.norm, reorderMeta=list(Condition=1), minMembers=minmembers)

if(diff_method=="DBA_DESEQ2"){
dba.ob.cont <- dba.analyze(dba.ob.cont, bBlacklist=F, bGreylist=F, method=DBA_DESEQ2)
#dba.ob.cont <- dba.analyze(dba.ob.cont, bBlacklist=F, bGreylist=F, method=DBA_DESEQ2)
} else if(diff_method=="DBA_EDGER") {
  dba.ob.cont <- dba.analyze(dba.ob.cont, bBlacklist=F, bGreylist=F, method=DBA_EDGER)
} else if(diff_method=="DBA_ALL_METHODS"){
  dba.ob.cont <- dba.analyze(dba.ob.cont, bBlacklist=F, bGreylist=F, method=DBA_ALL_METHODS)
}

#+ Fig3, fig.cap = "Fig 3: Correlation heatmap, using only significantly differentially bound sites" , fig.align = "center"
plot(dba.ob.cont,contrast=1)

dba.ob.diff <- dba.report(dba.ob.cont)


#+ Fig4, fig.cap = "Fig 4: Ven diagram of Gain vs Loss differentially bound sites" , fig.align = "center"
dba.plotVenn(dba.ob.cont,contrast=1,bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE)

#+ Fig5, fig.cap = "Fig 5: MA plot of Resistant-Responsive contrast. Sites identified as significantly differentially bound shown in red" , fig.align = "center"
dba.plotMA(dba.ob.cont)

#+ Fig6, fig.cap = "Fig 6: Volcano plot of Resistant-Responsive contrast. Sites identified as significantly differentially bound shown in red" , fig.align = "center"
dba.plotVolcano(dba.ob.cont)

write.table(dba.ob.diff, file=out_path, sep="\t", col.names=T, row.names=F, quote=F)

#' Below table displays first 6 rows of the differential analysis results using `r diff_method`
head(dba.ob.diff)
print("complted differential binding analysis")

