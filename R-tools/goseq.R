# Perform a GO term analysis with package goseq
# Maxime Caron Nov 2011
# Borrowed some argument code from Mathieu Bourgey
# usage : Rscript goseq.R -d path_dge_results_file -c comma_delimited_columns -g organism -k known_reference -o output_dir
# modified  by Mathieu Bourgey Feb 2013
# modified  by Mathieu Bourgey July 2013


usage=function(errM) {
        cat("\nUsage : Rscript goseq.R [option] <Value>\n")
        cat("       -d        : file containing the dge results\n")
        cat("       -c        : columns to use for input (first column = geneSymbol list, second column= adjusted p-values\n")
        cat("       -f        : FDR threshold detection of pathway enrichement\n")
        cat("       -p        : p-value threshold for e feature to be included in the GO analysis\n")
	cat("       -k        : known genome gene name mapping file\n")
	cat("       -m        : maximum number of pathway return (default all)\n")
	cat("       -v        : DGE method 0: edger/deseq ; 1: cuffdiff (default 0)\n")
	cat("       -i        : gene Name ID (default geneSymbol)\n")
	cat("       -G        : Nonnative Go file path\n")
	cat("       -a        : Gene length file path\n")
	cat("       -o        : output directory\n")
        cat("       -h        : this help\n\n")
        stop(errM)
}


##################################
ARG = commandArgs(trailingOnly = T)


## defult arg values
file=""
columns=""
fdrThr=0.05
pvalThr=0.05
out_path=""
maxP= -1
method=0
gene_path=NULL
go_path=""
gene_type="geneSymbol"

## get arg variables
for (i in 1:length(ARG)) {
        if (ARG[i] == "-d") {
                file=ARG[i+1]
        } else if (ARG[i] == "-c") {
                columns=ARG[i+1]
        } else if (ARG[i] == "-f") {
                fdrThr=as.numeric(ARG[i+1])
	} else if (ARG[i] == "-p") {
                pvalThr=as.numeric(ARG[i+1])
	} else if (ARG[i] == "-o") {
                out_path=ARG[i+1]
	} else if (ARG[i] == "-k") {
                known_ref=ARG[i+1]
	} else if (ARG[i] == "-m") {
                maxP=as.numeric(ARG[i+1])
	} else if (ARG[i] == "-v") {
                method=as.numeric(ARG[i+1])
	} else if (ARG[i] == "-G") {
                go_path=ARG[i+1]
	} else if (ARG[i] == "-a") {
                gene_path=ARG[i+1]
	} else if (ARG[i] == "-i") {
                gene_type=ARG[i+1]
	} else if (ARG[i] == "-h") {
                usage("")
        }
}
## check arg consitency
if (!(file.exists(file))) {
	stop("Input file not found") 
}
if (out_path == "") {
	stop("Output directory not found")
}


library('goseq')


print(paste("Using Non-native Gene Identifier",gene_path,"and category test",go_path))


set.seed(123456789)

d1<-read.table(file, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="")
toUse<-unlist(strsplit(columns,","))
selecT<-c(as.numeric(toUse[1]),as.numeric(toUse[2]))
d2<-d1[,selecT]

if(method == 1) {
kgX<-read.table(known_ref, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="")
tmp<-merge(d2,kgX,by.x=1,by.y=1)
d2<-tmp[,c(3,2)]
}
#data_externe[order(data_externe[,5]),]
d2<-d2[order(d2[,2]),]

head(d2)
is.significant<-function(x,pv=pvalThr) ifelse(x <= pv,1,0)
if(sum(is.significant(d2[,2])==1) == 0) {
	print("No significant adjusted p-values found")
	write.table(paste("Enriched category",paste("FDR <",as.character(fdrThr),"filtered p-value"),"GOID","Term","Ontology","Definition","Synonym", sep="\t"), out_path, append=F, row.names=F, col.names=F, quote=F)
    q("no",0)
}
d3<-cbind(d2[,1], is.significant(d2[,2]))
de<-subset(d3,d3[,2]==1)
gene.vector = as.integer(unique(d3[,1]) %in% unique(de[,1]))
names(gene.vector) = unique(d3[,1])
goTable=read.table(go_path,header=F, sep="\t", quote="", stringsAsFactors=F, comment.char="")
geneLenPre=read.table(gene_path,header=F, sep="\t", quote="", stringsAsFactors=F, comment.char="")
geneTable=unfactor(geneLenPre[,2])
names(geneTable)=unfactor(geneLenPre[,1])
pwf = nullp(gene.vector[names(gene.vector) %in% names(geneTable)], bias.data=geneTable[names(geneTable) %in% names(gene.vector)])
GO.wall =  goseq(pwf,gene2cat = goTable)


head(GO.wall)
enriched.GO = cbind(GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") <= fdrThr], GO.wall$over_represented_pvalue[p.adjust(GO.wall$over_represented_pvalue, method = "BH") <= fdrThr])
head(enriched.GO)
library(GO.db)
if(dim(enriched.GO)[1] == 0) {
	print(paste("No significant (FDR <",as.character(fdrThr),") enriched categories"))
	write.table(paste("Enriched category",paste("FDR <",as.character(fdrThr),"filtered p-value"),"GOID","Term","Ontology","Definition","Synonym", sep="\t"), out_path, append=F, row.names=F, col.names=F, quote=F)
	q("no",0)
} else {
if (maxP == -1) maxP=dim(enriched.GO)[1]
}
#write.table("Gene Ontology Analysis", out_path, row.names=F, col.names=F, quote=F)
write.table(paste("Enriched category",paste("FDR <",as.character(fdrThr),"filtered p-value"),"GOID","Term","Ontology","Definition","Synonym", sep="\t"), out_path, append=F, row.names=F, col.names=F, quote=F)
for (i in 1:maxP) {
f<-GOTERM[[enriched.GO[i,1]]]
if (!is.null(f)){
	write.table(paste(i, enriched.GO[i,2], GOID(f), Term(f), Ontology(f), Definition(f), Synonym(f)[i], sep="\t"), out_path, append=T, row.names=F, col.names=F, quote=F)
}
}

#getgo(d3[1:1000,15], "mm9", "geneSymbol", fetch.cats=c("GO:CC","GO:BP","GO:MF"))
