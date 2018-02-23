#!/usr/bin/env Rscript

##################################################################
##                                                              ##
##                           RobusTAD.R                         ##  
##       Rola Dali, Guillaume Bourque, Mathieu Blanchette       ##
##                           Feb 2018                           ##
##                             v1.0                             ## 
##                                                              ##
##################################################################

## The following script produces TAD Boundary scores for each bin on a chromosome.
## Input: interaction frequency matrix for a chromosome
## Output: 2 files: I- file with TAD Boundary scores: Right Boundary scores, Left Boundary scores and Final combined score (max(R, L)).
##                 II- file with TAD boundary calls identified by looking for local maxima above the set threshold


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ settings:


options(scipen=999)
library("optparse")

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ functions:


leftRightBoundaryScores <- function(data, minW, maxW, minRatio) {
  ## function that calculates left and right boundary scores for each bin
  s = dim(data)[1]
  datanorm <- matrix(0,dim(data)[1],dim(data)[2])
  
  if(status == "raw"){
    ## normalize if data is raw, skip otherwise
    print("... ... Normalizing Data")
    sumrow <- matrix(0,dim(data)[1],1)
    for (i in 1:s) {
      sumrow[i,1] = sum(data[i,])
    }
    ssum = sum(sumrow)
    
    for (i in 1:s) {
      for (j in 1:s) {
        datanorm[i,j]=data[i,j]/(sumrow[i,1]*sumrow[j,1])*ssum;
      }
    }
  } else {
    print("... ... skipping normalization step!")
    datanorm <- data
  }
  
  print("... ... Calculating Boundary Scores")
  leftBoundaryScore=matrix(1,s);
  rightBoundaryScore=matrix(1,s);
  
  for (i in (maxW+2):(s-maxW-2)) {
    rightScore=0;
    leftScore =0;
    leftBoundaryScore[i+1,1]=-99999;
    rightBoundaryScore[i+1,1]=-99999;
    
    diag=matrix(NA,maxW,maxW);
    for (d in 1:maxW) {
      maxiPossible = (d^3 + 3 * d^2 + 2* d)/6;
      for (j in 1:d) {
        diag[d,j] = datanorm[i-j+1,i-j+1+d];
      }
      
      for (del in 1:d) {
        diff = sum(datanorm[i-d,i-d+del]/(diag[del,]+0.001)>minRatio, na.rm=T) - sum(datanorm[i-d,i-d+del]/(diag[del,]+0.001)<1.0/minRatio, na.rm=T)
        rightScore = rightScore + diff;
      }
      if (d>=minW) {
        z = rightScore/sqrt(maxiPossible);
        # note: the score calculated pertains to the presence of a boundary on the right end of bin i. However it is more convenient to report scores relating to the left end of bins, so we associate this score to bin i+1 instead of i.
        if (z>rightBoundaryScore[i+1,1]) {
          rightBoundaryScore[i+1,1]=z;
        }
      }
      
      for (del in 1:d) {
        diff = sum(datanorm[i+1+d-del,i+1+d]/(diag[del,]+0.001)>minRatio, na.rm=T) - sum(datanorm[i+1+d-del,i+1+d]/(diag[del,]+0.001)<1.0/minRatio, na.rm=T)
        leftScore = leftScore + diff;
      }
      if (d>=minW) {
        z = leftScore/sqrt(maxiPossible);
        # note: the score calculated pertains to the presence of a boundary on the right end of bin i. However it is more convenient to report scores relating to the left end of bins, so we associate this score to bin i+1 instead of i.
        if (z>leftBoundaryScore[i+1,1]) leftBoundaryScore[i+1,1]=z;
      }
    }
  }
  return(cbind(leftBoundaryScore,rightBoundaryScore));
}


localpeaks <- function(score){
  ## locates peaks across the TAD score profile
  index <- c()
  for (i in 3:(length(score)-2)){
    if(score[i] > score[i-2] && score[i] > score[i-1] && score[i] > score[i+1] && score[i] > score[i+2]){
      index <- c(index, i)
    } 
  }
  return(index)
}



TADBoundcalls <- function(boundaries, threshold){
  ## calls TAD Boundaries using a threshold T by locating the maxima that are above the threshold
  BoundIndexR <- localpeaks(boundaries[,3])
  BoundIndexL <- localpeaks(boundaries[,2])
  BoundIndex <- unique(sort(c(BoundIndexR, BoundIndexL)))
  Boundmaxima <- boundaries[BoundIndex,]
  TADcalls <- Boundmaxima[Boundmaxima[,4] >= threshold,]
  return(TADcalls)
}


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ main:


## parse command line arguments:

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Interaction Frequency Matrix. Must be a square matrix: number of columns = number of rows"),
  make_option(c("-H", "--header"), action = "store_true", default=TRUE, 
              help="include -H if input contains a header/column names"),
  make_option(c("-n", "--norm"), action = "store", default ="raw", 
              help="indicates if IF matrix is raw or normalized [default = %default]; [options: {raw, norm}]"),
  make_option(c("-o", "--outDir"), type="character", default=".", 
              help="output directory name"),
  make_option(c("-b", "--binsize"), type="integer", default=50, 
              help="binsize or resolution used in Hi-C analysis in kb [default = %default]"),
  make_option(c("-r", "--minRatio"), type="double", default=1.5, 
              help="minimum ratio of Within to Across IF values to contribute to boundary score calculation [default = %default]"),
  make_option(c("-w", "--minWin"), type="integer", default=250, 
              help="minimum window around the bin used to calculate the TAD score in kb [default = %default]"),
  make_option(c("-W", "--maxWin"), type="integer", default=500, 
              help="maximum window around the bin used to calculate the TAD score in kb [default = %default]"),
  make_option(c("-T", "--threshold"), type="double", default=0.2, 
              help="data percentile of TAD scores used to calculate threshold in order to call significant TAD boundaries. [default = %default]; [options: 0-1];\n\t\tthe lower the threshold, the more stringent the TAD calls")
); 

description = "\n\n==============================================================================================================================================================\n\nRobusTAD calculates TAD Boundary scores for each bin on a chromosome.\n\tInput: interaction frequency matrix for a chromosome.\n\tOutput: 2 files: \n\t\t\tI- file with TAD boundary scores (BoundaryScores_*): contains Right Boundary scores, Left Boundary scores and Final combined score (max(R, L)).\n\t\t\tII- file with TAD boundary calls (TADBoundaryCalls_*) identified by looking for local maxima above the set threshod.\n\nUsage: RobusTAD.R -i InputMatrix [options]\n\n==============================================================================================================================================================\n\n"

epilogue = "\n\n==============================================================================================================================================================\n\nRobusTAD is available under a GPL liscence and comes with no warranties @ https://github.com/rdali/RobusTAD\n\n==============================================================================================================================================================\n\n"

opt_parser = OptionParser(option_list=option_list, description = description, epilogue=epilogue);
opt = parse_args(opt_parser);


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Error! missing arguments... Interaction frequency matrix is required!\n", call.=FALSE)
}


mymatrix <- opt$input
minRatio <- as.numeric(opt$minRatio)
binSize <- as.numeric(opt$binsize)
minWin <- as.numeric(opt$minWin)
maxWin <- as.numeric(opt$maxWin)
perc <- as.numeric(opt$threshold)
outDir <- opt$outDir
status <- opt$norm
header <- opt$header


if(!file.exists(mymatrix)){
  stop ("Error! provided matrix file does not exist")
}

##create outDir if it doesnt exist:
if (!file.exists(outDir)) {
  dir.create(outDir)
}


if(status != "raw" & status != "norm" ){
  stop ("Error! -n must be set to 'raw' or 'norm'")
}

if(perc < 0 | perc > 1){
  stop ("Error! --threshold should be set to a value between 0 and 1!")
}

minW <- as.integer(minWin/binSize)
maxW <- as.integer(maxWin/binSize)


## extract base name of matrix to use in output files:
prefix <- basename(mymatrix)
prefix <- substr(prefix, 1, nchar(prefix)-4)


## read in matrix
if(header){
  data <- as.matrix(read.table(mymatrix, header = TRUE, sep="\t", row.names = 1, as.is=TRUE))
} else {
  data <- as.matrix(read.table(mymatrix, header = FALSE, sep="\t", as.is=TRUE))
}

## check that input file is a square matrix:
if(nrow(data) != ncol(data)){
  stop ("Error! provided matrix is not a square matrix!")
}

boundaries = leftRightBoundaryScores(data, minW, maxW, minRatio)
rownames(boundaries) <- rownames(data)
colnames(boundaries) <- c("LeftBoundaryScore", "RightBoundaryScore")


## scale data to produce z-scores:
boundaries <- scale(boundaries, center = T, scale = T)

## calculate combined TAD score:
TADscore <- apply(boundaries, 1, max)
boundaries <- data.frame(boundaries, TADscore)

coordinates <- rownames(boundaries)
boundaries <- data.frame(boundaries, coordinates)
boundaries <- boundaries[,c(4,1,2,3)]


## call TAD boundaries:
print("... ... Calling Significant TAD Boundaries")
threshold <- quantile(boundaries$TADscore, na.rm = T, probs = (1-perc))[[1]]
TADcalls <- TADBoundcalls(boundaries, threshold)


## Save files:
write.table(boundaries, file.path(outDir, paste0("BoundaryScores_", prefix, "_binSize",binSize ,"_minW", minWin, "_maxW", maxWin,"_minRatio", minRatio,".txt")), quote = F, row.names=F)
write.table(TADcalls, file.path(outDir, paste0("TADBoundaryCalls_", prefix, "_binSize",binSize ,"_minW", minWin, "_maxW", maxWin,"_minRatio", minRatio,"_threshold",perc, ".txt")), quote = F, row.names=F)

print("... ... RobusTAD run is complete!")


