# Runs PopSV pipeline
# Written by Jean Monlong - Aug 2018

library(PopSV)
library(methods)

initfilename <- function(sample_bam_file, popsv_config_file){
  message('Create configuration file and folders')
  bam.files = read.table(sample_bam_file, as.is=TRUE, header=TRUE)
  files.df = init.filenames(bam.files)
  save(files.df, file=popsv_config_file)
  message('Done')
}

prepbins <- function(bin_file_input, bin_size, nb_chunks, bin_file){
  message('Prepare bin file')
  if(file.exists(bin_file_input)){
    message('File exists. Reading ', bin_file_input)
    bins.df = read.table(bin_file, as.is=TRUE, header=TRUE, sep='\t')
  } else {
    message('Binning hg19 in windows of ', bin_size, 'bp')
    bins.df = fragment.genome.hg19(bin_size)
  }
  if(!('GCcontent' %in% colnames(bins.df))){
    message('Retrieve GC content (hg19)')
    bins.df = getGC.hg19(bins.df)
  }
  message('Prepare chunks of bins')
  sm.chunk.size = ceiling(nrow(bins.df)/nb_chunks)
  bins.df = chunk.bin(bins.df, bg.chunk.size=5e5, sm.chunk.size=sm.chunk.size)
  write.table(bins.df, file=bin_file)
  message('Done')
}

countsample <- function(samp_name, popsv_config_file, bin_file){
  message('Count reads for ', samp_name)
  samp_name = make.names(samp_name)
  load(bins_file)
  load(popsv_config_file)
  file.i = which(files.df$sample == samp_name)
  bam.f = files.df$bam[file.i]
  bb.o = bin.bam(bam.f, bins.df, files.df$bc[file.i])
  message('Correct for GC bias')
  correct.GC(files.df$bc.gz[file.i], bins.df, files.df$bc.gc[file.i])
  message('Done')
}

preprefs <- function(popsv_config_file, ref_samps_file, ref_file, cont_sample_file, bins_file, graph_out, max_nb_refs){
  load(bins_file)
  load(popsv_config_file)
  message('Read the list of reference sample names')
  ref.samps = make.names(scan(ref_samps_file, ''))
  files.df = files.df[which(files.df$sample %in% ref.samps),]
  message('Merge read count for reference samples')
  pdf(graph_out)
  qc.o = qc.samples(files.df, bins.df, ref_file, nb.ref.samples=max_nb_refs)
  ## nb.cores=6 HOW TO SPECIFIY?
  dev.off()
  write(qc.o$cont.sample, file=cont_sample_file)
  message('Done')
}

normrefs <- function(ref_file, bin_file, cont_sample_file, chunk, res_file){
  load(bin_file)
  chunks = sort(unique(as.character(bins.df$sm.chunk)))
  chunk = chunks[chunk]
  bins.df.chunk = bins.df[which(bins.df$sm.chunk==chunk),]
  bg.chunk = head(bins.df.chunk$bg.chunk, 1)
  bc.df = read.bedix(ref_file, bins.df[which(bins.df$bg.chunk==bg.chunk),])
  cont.sample = scan(cont_sample_file, '')
  tn.norm(bc.df, cont.sample, bins=bins.df.chunk$bin)
}

mergeoutrefs <- function(popsv_config_file, norm_ref_prefix, nb_chunks, ref_prefix){

}

callsample <- function(samp_name, popsv_config_file, bin_file, cnv_file){
  samp_name = make.names(samp_name)
}


## Read arguments and call function
ARG = commandArgs(trailingOnly = T)
if(ARG[1] == 'initfilename'){
  initfilename(ARG[2], ARG[3])
} else if(ARG[1] == 'prepbins'){
  prepbins(ARG[2], ARG[3], ARG[4])
} else if(ARG[1] == 'countsample'){
  countsample(ARG[2], ARG[3], ARG[4])
} else if(ARG[1] == 'preprefs'){
  preprefs(ARG[2], ARG[3], ARG[4], ARG[5], ARG[6], ARG[7], ARG[8])
} else if(ARG[1] == 'normrefs'){
  normrefs(ARG[2], ARG[3], ARG[4], ARG[5], ARG[6])
} else if(ARG[1] == 'mergeoutrefs'){
  mergeoutrefs(ARG[2], ARG[3], ARG[4], ARG[5])
} else if(ARG[1] == 'callsample'){
  callsample(ARG[2], ARG[3], ARG[4], ARG[5])
} else {
  stop('Unknown step: ', ARG[1])
}


## TODO/QUESTIONS
# hg19 vs GRCh38
# multi-core jobs specification
# make.names and pb with sample name conversion
