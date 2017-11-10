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

##Extract the information from the seqz file
test <- sequenza.extract(data.file, parallel=3)

##Inference of cellularity and ploidy
CP.example <- sequenza.fit(test)

##Results of model fitting
sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = sample_name, out.dir=outputDir)

##detailed plot Confidence intervals, confidence region and point estimate
cint <- get.ci(CP.example)
cellularity <- cint$max.cellularity
ploidy <- cint$max.ploidy
avg.depth.ratio <- mean(test$gc$adj[, 2])

write(paste(sample_name,as.character(cellularity),as.character(ploidy),sep="\t"),paste(outputDir,paste(sample_name,"ploidy_celularity.tsv",sep="_"),sep="/"))

pdf(paste(outputDir,paste(sample_name,"ploidy_celularity.pdf",sep="_"),sep="/"))
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

