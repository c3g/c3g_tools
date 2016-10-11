## sequenza is not compatible with b37 natively
## we need to pre-select chromosme 1:22 + Y + X
## Mathieu bourgey - mathieu.bourgey@mcgill.ca
## 2015/12/02
## usage:
##   Rscript RunSequenza_analysis.R <FILE_sqeqz.hg19.gz> <OUTPUT_FOLDER> <OUTPUT_BASE_NAME>

##get args
args=commandArgs(TRUE)
data.file=args[1]
outputDir=args[2]
sample_name=args[3]

##load sequenza
library("sequenza")

##load result of sequenza-utils.py
seqz.data <- read.seqz(data.file)

##GC normalization
gc.stats <- gc.norm(x = seqz.data$depth.ratio, gc = seqz.data$GC.percent)
gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
seqz.data$adjusted.ratio <- seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)]

pdf(paste(outputDir,paste(sample_name,"GCnomr.pdf",sep="."),sep="/"))
par(mfrow = c(1,2), cex = 1, las = 1, bty = 'l')
matplot(gc.stats$gc.values, gc.stats$raw, type = 'b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2), xlab = 'GC content (%)', ylab = 'Uncorrected depth ratio')
legend('topright', legend = colnames(gc.stats$raw), pch = c(1, 19, 1))
hist2(seqz.data$depth.ratio, seqz.data$adjusted.ratio, breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2), xlab = 'Uncorrected depth ratio', ylab = 'GC-adjusted depth ratio')
dev.off()

##Extract the information from the seqz file
test <- sequenza.extract(data.file)


##Inference of cellularity and ploidy
CP.example <- sequenza.fit(test)

##Results of model fitting
sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = sample_name, out.dir=outputDir)

##detailed plot Confidence intervals, confidence region and point estimate
cint <- get.ci(CP.example)
cellularity <- cint$max.cellularity
ploidy <- cint$max.ploidy
avg.depth.ratio <- mean(test$gc$adj[, 2])

write(paste(sample_name,as.character(cellularity),as.character(ploidy),sep="\t"),paste(outputDir,paste(sample_name,"Ploidy_celularity.tsv",sep="."),sep="/"))

pdf(paste(outputDir,paste(sample_name,"Ploidy_celularity.pdf",sep="."),sep="/"))
par(mfrow = c(2,2))
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE,col='red')
plot(cint$values.cellularity, ylab = "Cellularity",xlab = "posterior probability", type = "n")
select <- cint$confint.cellularity[1] <= cint$values.cellularity[,2] & cint$values.cellularity[,2] <= cint$confint.cellularity[2]
polygon(y = c(cint$confint.cellularity[1], cint$values.cellularity[select, 2], cint$confint.cellularity[1]) ,x = c(0, cint$values.cellularity[select, 1], 0), col='red', border=NA)
lines(cint$values.cellularity)
abline(h = cint$max.cellularity, lty = 2, lwd = 0.5)
plot(cint$values.ploidy, xlab = "Ploidy", ylab = "posterior probability", type = "n")
select <- cint$confint.ploidy[1] <= cint$values.ploidy[,1] & cint$values.ploidy[,1] <= cint$confint.ploidy[2]
polygon(x = c(cint$confint.ploidy[1], cint$values.ploidy[select, 1], cint$confint.ploidy[1]), y = c(0, cint$values.ploidy[select, 2], 0), col='red', border=NA)
lines(cint$values.ploidy)
abline(v = cint$max.ploidy, lty = 2, lwd = 0.5)
mtext(paste("Estimated Cellularity:",as.character(cellularity),"and Ploidy:",as.character(ploidy)), side = 3, line=-2, outer = TRUE)
dev.off()

##Detect variant alleles (mutations)
mut.tab = na.exclude(do.call(rbind, test$mutations))
mut.alleles = mufreq.bayes(mufreq = mut.tab$F, depth.ratio = mut.tab$adjusted.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)

somatic.mut = cbind(mut.tab[,c("chromosome","position","F","adjusted.ratio", "mutation")],mut.alleles)
write.table(somatic.mut,paste(outputDir,paste(sample_name,"Somatic_mutations.tsv",sep="."),sep="/"),quote=F,row.names=F,sep="\t")

##Detect copy number variations
seg.tab = na.exclude(do.call(rbind, test$segments))
cn.alleles = baf.bayes(Bf = seg.tab$Bf, depth.ratio = seg.tab$depth.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)
seg.tab <- cbind(seg.tab, cn.alleles)
write.table(seg.tab,paste(outputDir,paste(sample_name,"Somatic_CNsegments.tsv",sep="."),sep="/"),quote=F,row.names=F,sep="\t")

