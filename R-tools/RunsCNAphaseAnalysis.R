##get args
args=commandArgs(TRUE)
file_name=args[1]
normal_name=args[2]
tumor_name=args[3]

##load sCNAphase
library("sCNAphase") 

chroms=c(1:22)

##run sCNAphase and generate purity and plodiy estimates
#   This will generate a R dat file in the current directory called res.{anaName}.phased.chr.W.dat, based on which sCNAphase can then format the estimation to the segmentation file, the d.SKY plot and a vcf file.

inferCNA(file_name, normal_name, tumor_name, chroms, allelicMapability=T)

##Generates a *.csv file. Each row corresponds
#   to a particular sCNAs with chrID, start, end, copy number, allelic copy number.

genSegFile(anaList= file_name, outdir = "results")

##Generates the d.SKY plot into a pdf file.

produceDSKY(anaList= file_name, outDir = "results")
