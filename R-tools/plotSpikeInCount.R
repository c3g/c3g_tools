#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
# args[1] filter blast file
# args[2] total number of reads sampled

file_name=args[1]
data=read.delim(file_name, as.is=TRUE) #check that not empty?
tags=data[,"Tag"]
total_tags=nrow(data)
alignments=data[,"Alignment"]
identity=data[,"Identity"]
mismatches=data[,"Mismatches"]
samples_size=as.integer(args[2])

tags.table.count=table(tags)
tags.table.count=tags.table.count[order(tags.table.count, decreasing=T)]
tags.table=tags.table.count/samples_size*100
tag_color=c("lightblue",rainbow(length(tags.table)-1, alpha = 1, v=0.8))
alignments.table=table(tags,alignments)/total_tags*100
alignments.table=alignments.table[names(tags.table), ,drop=F] #to avoid losing one dimension if only 1 column
identity.table=table(tags,identity)/total_tags*100
identity.table=identity.table[names(tags.table), , drop=F]
mismatches.table=table(tags,mismatches)/total_tags*100
mismatches.table=mismatches.table[names(tags.table), ,drop=F]

# IC at 95% for binomial distribution
IC95_bn<-function(p, n){
	d=sqrt(p*(1-p)/n)*1.96
	ic=c(p-d, p+d)
	return(ic)}

###########################
#plot blast results
###########################
pdf(paste0(file_name,".pdf"))
if(nrow(data)>0) {
	tags.ic95.count=sapply(tags.table.count, function(x) {IC95_bn(x/samples_size, samples_size)})*samples_size
	tags.ic95=sapply(tags.table, function(x) {IC95_bn(x/100, samples_size)*100})
	m <- matrix(c(1,1,2,3,4,5),nrow = 3,ncol = 2, byrow = TRUE)
	layout(mat = m,heights = c(0.1,0.45,0.45))
	#legend
	par(mar = c(0,0,2,0))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="", main=basename(file_name),cex.main=1)
	legend(x = "top",inset = 0, legend = names(tags.table), fill=tag_color, border=tag_color, horiz = TRUE, title=paste0(samples_size, " sampled reads"))
	par(mar=c(4, 4, 4, 2) + 0.1)
	#tag count
	if (length(tags.table)>5){
		tags.bp<-barplot(tags.table.count,main="Tags", col=tag_color,border=tag_color, ylab="Read count", names.arg=tags.table.count, ylim=c(0, max(tags.table.count,tags.ic95.count[2,])), las=3)
	} else {
		tags.bp<-barplot(tags.table.count,main="Tags", col=tag_color,border=tag_color, ylab="Read count", names.arg=paste(tags.table.count, names(tags.table),sep="\n"), ylim=c(0, max(tags.table.count,tags.ic95.count[2,])))
	}
	suppressWarnings(arrows(tags.bp, tags.ic95.count[1,], tags.bp, tags.ic95.count[2,], angle=90, code=3, length=0.1))
	#tag proportion
	if (length(tags.table)>5){
		tags.bp<-barplot(tags.table,main="Tags", col=tag_color,border=tag_color, ylab="% of reads", names.arg=paste0(signif(tags.table,3),"%"), ylim=c(0, max(tags.table,tags.ic95[2,])), las=2)
	} else {
		tags.bp<-barplot(tags.table,main="Tags", col=tag_color,border=tag_color, ylab="% of reads", names.arg=paste(paste0(signif(tags.table,3),"%"), names(tags.table),sep="\n"), ylim=c(0, max(tags.table,tags.ic95[2,])))
	}
	suppressWarnings(arrows(tags.bp, tags.ic95[1,], tags.bp, tags.ic95[2,], angle=90, code=3, length=0.1))
	#aligment length
	barplot(alignments.table, main="Alignment length", col=tag_color, border=tag_color, ylab="% of alignments")
	#mismatches
	barplot(mismatches.table, main="Mismatches", col=tag_color, border=tag_color, ylab="% of alignments")
	par(mar=c(5, 4, 4, 2) + 0.1)
} else {
	m <- matrix(c(1,1,2,3,4,5),nrow = 3,ncol = 2, byrow = TRUE)
	layout(mat = m,heights = c(0.1,0.45,0.45))
	par(mar = c(0,0,2,0))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="", main=basename(file_name),cex.main=1)
	text(1,1,labels=c("No tag detected"))
	for (i in 1:4){
		plot(1, type = "n", axes=FALSE, xlab="", ylab="", main="",cex.main=0.7)
	}
}
out <- dev.off()

###########################
#print info summary
###########################
tags.summary=data.frame(matrix(NA, ncol=4, nrow=length(tags.table)))
colnames(tags.summary)=c("Tag", "Count", "Estimate_%" ,"IC95")
if(nrow(data)>0) {
	tags.summary[,"Tag"]=names(tags.table)
	tags.summary[,"Count"]=tags.table.count
	tags.summary[,"Estimate_%"]=tags.table
	tags.summary[,"IC95"]=paste0("[",signif(tags.ic95[1,],3),",", signif(tags.ic95[2,],3),"]")
}
write.table(tags.summary, file=paste0(file_name,".summary"), quote=F, sep="\t", col.names=T, row.names=F)

###########################
# print tag distribution
###########################
ggplot2 <- require(ggplot2)
if (ggplot2){
	out_pdf=paste0(file_name,"tagDistribution.pdf")

	read_length=NA
	biostrings=require("Biostrings")
	fasta_file=gsub("blast.*","fasta",file_name)
	if(biostrings & file.exists(fasta_file)){
		fasta=readDNAStringSet(fasta_file)
		if(length(unique(width(fasta)))==1){
			read_length=unique(width(fasta))
		} else{
			warning("Non uniform read length, tag distribution might be incorrect.")
		}
	}
	if(is.na(read_length)){
		read_length=max(data[,"Query_end"]) #assume all reads have the same length!
		warning(paste("Setting read length to", read_length, "bp."))
	}

	pdf_out=paste0(file_name,".tagDistribution.pdf")
	data_query=apply(data[,c("Tag","Query_start","Query_end")], 1, function(x) {data.frame(Tag=x[[1]],Pos=seq(as.integer(x[[2]]),as.integer(x[[3]])))}) #create a list of 1's for each position with a tag
	data_query=do.call("rbind", data_query)
	p <- ggplot(data_query, aes(Pos, fill=Tag)) +  geom_histogram(binwidth = 1, show.legend=FALSE) + xlim(0, read_length) + facet_grid(Tag ~ ., scales="free_y") + scale_fill_manual(values=tag_color)
	# ggsave(out_pdf)
	p <- p + ggtitle("Distribution of tags along the reads") + ylab("Tag count") + xlab("Position")
	p
	ggsave(pdf_out)
}



print("PDF and summary files successfully generated")