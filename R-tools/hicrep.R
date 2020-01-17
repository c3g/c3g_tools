# Evaluates reproducibility of Hi-C data using hicrep
# Inputs: two hic-interaction matrices of same chromosme (from two different samples)
# Written by Pubudu Nawarathna - Nov 2019
# Usage : Rscript hicrep.R -s1 sample1_name -s2 sample2_name -f1 file_path1 -f2 file_path2 -c chromosome -o output_file_path -sm h/smooth_value -r resolution/bin size 
# -b bundary -w save_weights -cor save_correlation_matirx

library(data.table)
library(hicrep)

# Usage

usage = function(errM) {
  cat("\nUsage : Rscript hicrep.R [option] <Value>\n")
  cat("       -s1         : sample1\n")
  cat("       -s2         : sample2\n")
  cat("       -f1         : file1_path\n")
  cat("       -f2         : file2_path\n")
  cat("       -c          : chr\n")
  cat("       -o          : output file path\n")
  cat("       -sm         : smoothing variable h\n")
  cat("       -r          : bin size (resolution)\n")
  cat("       -b          : boundary\n")
  cat("       -w          : weights\n")
  cat("       -cor        : correlation\n")
  cat("       -h          : this help\n\n")
  
  stop(errM)
}

set.seed(123456789)

#The inputs for perform_hicrep function are simply two Hi-C matrices to be compared.
#The Hi-C matrices should have the dimension N×(3+N).
#The three additional initial columns are the chromosome name
#and mid-point coordinates of the two contact bins


perform_hicrep <- function(mat1=NULL, mat2=NULL, bin=1000000, smooth=5, boundary=5000000, restructure_data=T, down_sampling="TRUE" ){
  
  mat1[is.na(mat1)] <- 0;
  mat2[is.na(mat2)] <- 0;
  
  print(paste("matrix 1 size=", dim(mat1)))
  print(paste("matrix 2 size=", dim(mat2)))
  
  if(dim(mat1)[1]!=dim(mat2)[1]){
    if(dim(mat1)[1] < dim(mat2)[1]){
      mat2 <- mat2[1:dim(mat1)[1], 1:dim(mat1)[2]]
      print("matrix 2 size is adjusted")
    }
    else{
      mat1 <- mat1[1:dim(mat2)[1], 1:dim(mat2)[2]]
      print("matrix 1 size is adjusted")
    }
  }
  
  if(down_sampling=="TRUE"){
    
    #Down-sampling based on the smallest matrix
    #only the large matrix will be down-sampled
    print(paste("matrix 1 rows=",as.numeric(sum(mat1[,-c(1:3)]))))
    print(paste("matrix 2 rows=",as.numeric(sum(mat2[,-c(1:3)]))))
    
    if(as.numeric(sum(mat1[,-c(1:3)])) == as.numeric(sum(mat2[,-c(1:3)]))){
    } else if(as.numeric(sum(mat1[,-c(1:3)])) < as.numeric(sum(mat2[,-c(1:3)]))){
      size <- as.numeric(sum(mat1[,-c(1:3)])) 
      mat2_h <- depth.adj(mat2, size, as.numeric(bin), out = 0)
      print("matrix 2 down-sampled")
      mat1_h <- mat1
      
    }
    else{
      
      size <- as.numeric(sum(mat2[,-c(1:3)])) 
      mat1_h <-  depth.adj(mat1, size, as.numeric(bin), out = 0)
      print("matrix 1 down-sampled")
      mat2_h <- mat2
      
    }
    
  } else if(down_sampling=="FALSE"){
    
    #Doing nothing
    
  } else if(is.numeric(as.numeric(down_sampling))){
  
    mat2 <-  depth.adj(mat2, as.numeric(down_sampling), as.numeric(bin), out = 0)
    mat1 <-  depth.adj(mat1, as.numeric(down_sampling), as.numeric(bin), out = 0)
    print("both matrices were down-sampled")
    #Down-sampling based on the given value
    #both the samples will be down sampled
    
    mat1_h <- mat1
    mat2_h <- mat2
    
  }
  
  print(paste("sample 1 new rows=",as.numeric(sum(mat1[,-c(1:3)])) ))
  print(paste("sample 2 new rows=",as.numeric(sum(mat2[,-c(1:3)])) ))
  
  #if down_sampling==false no down-sampling will be performed
  
  ##if smooth==TRUE optimal smoothing parameter will search using htrain function
  if(smooth=="TRUE"){
  print("Calculaitng optimal smoothing paramter")
  smooth <- calculate_optimal_h(mat1_h, mat2_h, as.numeric(bin), as.numeric(boundary))
  print("smoothing value calculated")
  print(paste("Selected optimal smoothing value=",smooth))
  }
  Pre_HiC <- prep(mat1, mat2, as.numeric(bin), as.numeric(smooth), as.numeric(boundary))
  print("pcrep reproducibility score was calculated")
  SCC.out = get.scc(Pre_HiC, bin, boundary)
  return(list(SCC.out,smooth))
}


hicrep_analysis <- function(out_file=NULL, sample1=NULL,sample2=NULL, file1_path=NULL, file2_path=NULL, chr=NULL, bin=50000,smooth=5, boundary=500000,
                            down_sampling="TRUE", corr=FALSE, weights=FALSE){
  
 
  mat1 <- fread(file1_path, data.table=F, na.strings=c("",NA,"NULL"))
  mat2 <- fread(file2_path, data.table=F, na.strings=c("",NA,"NULL"))
  

  rscore <- perform_hicrep(mat1, mat2, bin=as.numeric(bin), boundary = as.numeric(boundary), down_sampling=down_sampling, smooth=smooth)
  
  
  write.table(paste(rscore[[1]][[3]], rscore[[1]][[4]] ,rscore[[2]][[1]], sep="\t"), file = out_file, row.names = F, col.names = F, quote = F)
  
  if(corr==TRUE) {
    write.table(as.data.frame(rscore[[1]][[1]]), file = paste0(out_file, ".corr"), row.names = F, col.names = T, quote = F)
  }
  if(weights==TRUE){
    write.table(as.data.frame(rscore[[1]][[2]]), file = paste0(out_file, ".weights"), row.names = F, col.names = T, quote = F)
  }

}

calculate_optimal_h <- function(mat1=NULL, mat2=NULL, bin=NULL, boundary=NULL){
  
  # A fraction (10%) of data are first randomly sampled, then the scc for the sampled data is computed
  # at a series of smoothing parameterts in the ascending order. The samllest h at which the increment
  # of scc is less than 0.01 is saved. This procedure is repeat 10 times, and the mode of the 10 h’s is
  # outputed as the estimated optimal neighborhood size.
  
   h_hat <- htrain(mat1, mat2, as.numeric(bin), as.numeric(boundary), 0:10)
   return(h_hat)
}


ARG = commandArgs(trailingOnly=T)

## default arg values
#no default args

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-s1") {
    sample1 = ARG[i+1]
  } else if (ARG[i] == "-s2") {
    sample2 = ARG[i+1]
  } else if (ARG[i] == "-o") {
    out_path = ARG[i+1]
  } else if (ARG[i] == "-c") {
    chr = ARG[i+1]
  } else if (ARG[i] == "-h") {
    usage("")
  } else if (ARG[i] == "-f1") {
    file1_path = ARG[i+1]
  } else if (ARG[i] == "-f2") {
    file2_path = ARG[i+1]
  } else if (ARG[i] == "-sm") {
    smooth = ARG[i+1]
  } else if (ARG[i] == "-r") {
    bin = ARG[i+1]
  } else if (ARG[i] == "-b") {
    boundary = ARG[i+1]
  } else if (ARG[i] == "-w") {
    weights = ARG[i+1]
  } else if (ARG[i] == "-cor") {
    corr = ARG[i+1]
  } else if (ARG[i] == "-d") {
    down_sampling = ARG[i+1]
  }
}

## check arg consitency
if (!(file.exists(file1_path))) {
  usage("Error : Sample1 file not found")
}
if (!(file.exists(file2_path))) {
  usage("Error : Sample2 file not found")
}
if (out_path == "") {
  usage("Error : Output directory not specified")
}
if (sample1 == "") {
  usage("Error : Name of the sample1 is not specified")
}
if (sample2 == "") {
  usage("Error : Name of the sample2 is not specified")
}
if (smooth == "") {
  usage("Error : smoothing variable h is not specified")
}
if (bin == "") {
  usage("Error : bin size is not specified")
}
if (boundary == "") {
  usage("Error : boundary size is not specified")
}
if (down_sampling == "") {
  usage("Error : down-sampling value is not specified")
}

print(paste0("sample1 name=",sample1))
print(paste0("sample1 file path=",file1_path))
print(paste0("sample2 name=",sample2))
print(paste0("sample2 file path=",file2_path))
print(paste("chromosome=",chr))
print(paste("smooth value calucaltion=",smooth))
print(paste("Down-sampling=",smooth))
print(paste("Boundary=",boundary))


hicrep_analysis(out_file=out_path, sample1=sample1, sample2=sample2, file1_path=file1_path, file2_path=file2_path, chr=chr, bin=bin,
                smooth=smooth, boundary=boundary, down_sampling=down_sampling, corr=as.logical(corr), weights=as.logical(weights))
