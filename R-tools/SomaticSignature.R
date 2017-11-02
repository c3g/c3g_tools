#load library works with module m mugqic_dev/R/3.2.2 mugqic_dev/R_Bioconductor/3.2.2_3.1
library(SomaticSignatures)

library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
library(Cairo)

#load somatic
files <- list.files("./",pattern=".vcf$",recursive=T,full.names=TRUE)
vranges <- lapply(files, function(v) readVcfAsVRanges(v,"hs37d5"))
vranges.cat <- do.call(c,vranges)
vranges.cat.tumor1=vranges.cat[grep(pattern="CUT",sampleNames(vranges.cat))]
vranges.cat.tumor2=vranges.cat[grep(pattern="-T",sampleNames(vranges.cat))]
vranges.cat.tumor= c(vranges.cat.tumor1,vranges.cat.tumor2)
sampleNames(vranges.cat.tumor)=factor(as.vector(sampleNames(vranges.cat.tumor)))

#Run the mutationContext function of SomaticSignatures. (takes time to run)
mc <- mutationContext(vranges.cat.tumor, BSgenome.Hsapiens.1000genomes.hs37d5)

#make a matrix that contains counts of mutations in each of the 96 possible combinations of mutations and contexts counting up the totals separately for each sample
mm <- motifMatrix(mc, group = "sampleNames", normalize=TRUE)

#find out how many signatures we have in the data run the command (didn't work for a single signature => 1:15)
gof_nmf <- assessNumberSignatures(mm, 2:10, nReplicates = 10)

#plot the results from the NMF
Cairo(file="plotNumberOfSignatures_NoOutlier.pdf", type="pdf", units="in", width=9, height=8, dpi=72)
x=plotNumberSignatures(gof_nmf)
print(x)
dev.off()

sigs_nmf_4 = identifySignatures(mm, 4, nmfDecomposition)

Cairo(file="plot_4Signatures_NoOutlier_class.pdf", type="pdf", units="in", width=10, height=8, dpi=72)
z=plotSignatures(sigs_nmf_4,normalize=TRUE, percent=FALSE) + ggtitle("Somatic Signatures: NMF - Barchart") + scale_fill_brewer(palette = "Set2")
show(z)
dev.off()

Cairo(file="PlotSampleContribution_4_Signatures.pdf", type="pdf", units="in", width=9, height=6, dpi=72)
y=plotSamples(sigs_nmf_4, normalize=TRUE) + scale_y_continuous(breaks=seq(0, 1, 0.2), expand = c(0,0))+ theme(axis.text.x = element_text(size=6))
show(y)
dev.off()

Cairo(file="PlotSampleContribution_4_heatmap.pdf", type="pdf", units="in", width=9, height=6, dpi=72)
pheatmap(samples(sigs_nmf_4),cluster_cols=F, clustering_distance_cols = "correlation")
dev.off()

#deconstructSigs
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

outSigTable=NULL
for (i in rownames(sigs.input)) {
	output.sigs = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.nature2013, sample.id = i)
        outSigTable=rbind(outSigTable,c(output.sigs$weights,output.sigs$unknown))
	#plot deconstructed signature
	Cairo(file=paste("PlotSample",i,"deconstructAlexandrov.pdf",sep="_"), type="pdf", units="in", width=9, height=6, dpi=72)
	plotSignatures(output.sigs)
	dev.off()
	# plot deconstructed signature pie
	Cairo(file=paste("PlotSample",i,"deconstructAlexandrov_pie.pdf",sep="_"), type="pdf", units="in", width=9, height=6, dpi=72)
	makePie(output.sigs)
	dev.off()
}
rownames(outSigTable)=rownames(sigs.input)
write.table(file="AllSample.signatureAlexamndrovfit.tsv", x=outSigTable,sep="\t",quote=F)


##fit with alexendrov signatures - decrepeted
library(entropy)
data(signatures21)
entroSig=matrix(rep(NA,21*7),ncol=7,nrow=21)
for (i in 1:7) {
    for (j in 1:21) {
        entroSig[j,i]=KL.empirical(mm[,i],signatures21[,j])
}
}
