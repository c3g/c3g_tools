# module load mugqic_dev/R_Bioconductor/3.2.2_3.1

# Get arguments
# Arg 1: file  with vcf uncompressed paths and tumor sample name
# Arg 2: output folder
args <- commandArgs(trailingOnly = TRUE)

# Load the analysis package
library(SomaticSignatures)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
library(Cairo)



files <- read.table(args[1],sep="\t")
vranges <- lapply(files$path, function(v) readVcfAsVRanges(v,"hs37d5"))
vranges.cat <- do.call(c,vranges)

###
vranges.cat.tumor=NULL
for (i in files$tumor) {
  vranges.cat.tumor=c(vranges.cat.tumor,vranges.cat[grep(pattern=i,sampleNames(vranges.cat))])
}

vranges.cat.tumor=vranges.cat.tumor[[1]]

sampleNames(vranges.cat.tumor)=factor(as.vector(sampleNames(vranges.cat.tumor)))
#Run the mutationContext function of SomaticSignatures. (takes time to run)
mc <- mutationContext(vranges.cat.tumor, BSgenome.Hsapiens.1000genomes.hs37d5)
#make a matrix that contains counts of mutations in each of the 96 possible combinations of mutations and contexts counting up the totals separately for each sample
mm <- motifMatrix(mc, group = "sampleNames", normalize=TRUE)

##
library(deconstructSigs)

sigs.input=as.data.frame(t(mm))
colnames(sigs.input)=c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G",
 "C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C",
 "T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A",
 "C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T",
 "T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G",
 "A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C",
 "G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A",
 "A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
 "G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G",
 "T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C",
 "C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A",
 "T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T",
 "C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G",
 "G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
#rownames(sigs.input)
AlexSigTable=NULL
for (i in rownames(sigs.input)) {
output.sigs = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.nature2013, sample.id = i)
AlexSigTable=rbind(AlexSigTable, output.sigs$weights)
#plot deconstructed signature
Cairo(file=file.path(arg[2],paste("PlotSample",i,"deconstructAlexandrov.pdf",sep="_")), type="pdf", units="in", width=9, height=6, dpi=72)
plotSignatures(output.sigs)
dev.off()
# plot deconstructed signature pie
Cairo(file=file.path(arg[2],paste("PlotSample",i,"deconstructAlexandrov_pie.pdf",sep="_")), type="pdf", units="in", width=9, height=6, dpi=72)
makePie(output.sigs)
dev.off()
}
write.table(AlexSigTable,file=file.path(arg[2],"Alexandrov_weigth.tsv"),quote=F,sep="\t")

