# Performs various filter of Axiam gene titan data
# Written by Mathieu Bourgey & Jean Molong - August 2016, April 2017
# Usage : Rscript filterAxiom.R -s path_design -c path_rawcountfile -o output_dir

# Usage functions

usageMain=function(errM) {
	cat("\nUsage : Rscript filterAxiom.R -s <Step>\n")
	cat("       genoqc    : filter sample based on Dish QC value (Best Practices Step 3) \n")
	cat("       sampleqc  : filter sample based on call rate value (Best Practices Step 5) \n")
	cat("       plateqc   : filter sample based on the Plate Pass Rat value (Best Practices Step 6) \n")
	cat("       merge     : merge apt output when the analysis is chunked \n")
	cat("       help      : this help\n\n")
	stop(errM)
}

usageGenoQC=function(errM) {
	cat("\nUsage : Rscript filterAxiom.R -s genoqc [option] <Value>\n")
	cat("       -c        : CEL files list to be filtered\n")
	cat("       -q        : apt-geno qc file\n")
	cat("       -d        : DQC lower value (default 0.82)\n")
	cat("       -o        : output file\n")
	cat("       -h        : this help\n\n")
	stop(errM)
}

usageSampleQC=function(errM) {
	cat("\nUsage : Rscript filterAxiom.R -s sampleqc [option] <Value>\n")
	cat("       -c        : CEL files list to be filtered\n")
	cat("       -q        : genotype summary file\n")
	cat("       -d        : minimal call rate value (default 97)\n")
	cat("       -o        : output file\n")
	cat("       -h        : this help\n\n")
	stop(errM)
}

usagePlateQC=function(errM) {
	cat("\nUsage : Rscript filterAxiom.R -s plateqc [option] <Value>\n")
	cat("       -q        : QC call rate file\n")
	cat("       -c        : CEL files list to be filtered\n")
	cat("       -i        : intermediate CEL files list after DQC step\n")
	cat("       -a        : minimal average call rate value of passing sample (default 98.5)\n")
	cat("       -d        : minimal pass plate value (default 95)\n")
	cat("       -p        : file matching CEL file and plates \n")
	cat("       -o        : output file\n")
	cat("       -h        : this help\n\n")
	stop(errM)
}

usagePlotSNP=function(errM) {
	cat("\nUsage : Rscript filterAxiom.R -s plotSNP [option] <Value>\n")
	cat("       -d        : number of SNP to plot per category (default 6)\n")
	cat("       -o        : output directory\n")
	cat("       -h        : this help\n\n")
	stop(errM)
}

usageMerge=function(errM) {
	cat("\nUsage : Rscript filterAxiom.R -s merge [option] <Value>\n")
	cat("       -n        : number of chunks (default 10)\n")
	cat("       -o        : output directory\n")
	cat("       -h        : this help\n\n")
	stop(errM)
}

##################################
main=function(){
	ARG = commandArgs(trailingOnly = T)

	## get arg variables
	if (length(ARG) > 1) {
		for (i in 1:length(ARG)) {
			if (ARG[i] == "-s") {
				if (ARG[i+1] == "genoqc") {
					genoQCmain(ARG)
					break
				} else if (ARG[i+1] == "sampleqc") {
					sampleQCmain(ARG)
					break
				} else if (ARG[i+1] == "plateqc") {
					plateQCmain(ARG)
					break
				} else if (ARG[i+1] == "plotSNP") {
					plotSNPmain(ARG)
					break
				} else if (ARG[i+1] == "merge") {
					mergeOutput(ARG)
					break
				} else {
					usageMain("Step not found !")
				}
			} else if (ARG[i] == "-h" || ARG[i] == "--help" ) {
				usageMain("")
			} else {
				usageMain("Step not found 2")
			}
		}
	} else {
		usageMain("MISSING ARGUMENT")
	}
}

#Best-practice step 3
genoQCmain=function(ARG) {
	#get specific args
	min_DQC=0.82
	out_file=""
	cel_file=""
	qc_file=""
	for (i in 1:length(ARG)) {
		if (ARG[i] == "-d") {
			min_DQC=as.numeric(ARG[i+1])
		} else if (ARG[i] == "-q") {
			qc_file=ARG[i+1]
		} else if (ARG[i] == "-c") {
			cel_file=ARG[i+1]
		} else if (ARG[i] == "-o") {
			out_file=ARG[i+1]
		} else if (ARG[i] == "-h") {
			usageGenoQC("")
		}
	}
	if (!(file.exists(cel_file))) {
		usageGenoQC("Error : CEL file list file not found")
	}
	if (!(file.exists(qc_file))) {
		usageGenoQC("Error : apt-geno qc file not found")
	}
	if (!(is.numeric(min_DQC))) {
		usageGenoQC("Error : threshold min DQC should be numeric [0-1]")
	}
	if (out_file == "" ) {
		usageGenoQC("Error : output CEL file list not specified")
	}
	celList=read.table(cel_file,header=T)
	qc_table=read.table(qc_file,header=T)
	sample_to_keep=qc_table$cel_files[qc_table$axiom_dishqc_DQC >= min_DQC]
	celList_filtered=data.frame(cel_files=celList$cel_files[basename(as.vector(celList$cel_files)) %in% sample_to_keep])
	write.table(celList_filtered,out_file,col.names=T,row.names=F,quote=F)
}

#Best-practice step 5
sampleQCmain=function(ARG) {
	#get specific args
	min_CR=97
	out_file=""
	cel_file=""
	cr_file=""
	for (i in 1:length(ARG)) {
		if (ARG[i] == "-d") {
			min_CR=as.numeric(ARG[i+1])
		} else if (ARG[i] == "-q") {
			cr_file=ARG[i+1]
		} else if (ARG[i] == "-c") {
			cel_file=ARG[i+1]
		} else if (ARG[i] == "-o") {
			out_file=ARG[i+1]
		} else if (ARG[i] == "-h") {
			usageSampleQC("")
		}
	}
	if (!(file.exists(cel_file))) {
		usageSampleQC("Error : CEL file list file not found")
	}
	if (!(file.exists(cr_file))) {
		usageSampleQC("Error : genotype call rate file not found")
	}
	if (!(is.numeric(min_CR))) {
		usageSampleQC("Error : threshold min sample QC should be numeric [0-100]")
	}
	if (out_file == "" ) {
		usageSampleQC("Error : output CEL file list not specified")
	}
	celList=read.table(cel_file,header=T)
	cr_table=read.table(cr_file,header=T)
	sample_to_keep=cr_table$cel_files[cr_table$call_rate >= min_CR]
	celList_filtered=data.frame(cel_files=celList$cel_files[basename(as.vector(celList$cel_files)) %in% sample_to_keep])
	write.table(celList_filtered,out_file,col.names=T,row.names=F,quote=F)
}

#Best-practice step 6
plateQCmain=function(ARG) {
	#get specific args
	min_CR=98.5
	min_PPR=95
	out_file=""
	cel_file=""
	cel_file2=""
	match_file=""
	cr_file=""
	dqc_file=""
	for (i in 1:length(ARG)) {
		if (ARG[i] == "-a") {
			min_CR=as.numeric(ARG[i+1])
		} else if (ARG[i] == "-q") {
			cr_file=ARG[i+1]
		} else if (ARG[i] == "-t") {
			dqc_file=ARG[i+1]
		} else if (ARG[i] == "-d") {
			min_PPR=as.numeric(ARG[i+1])
		} else if (ARG[i] == "-p") {
			match_file=ARG[i+1]
		} else if (ARG[i] == "-c") {
			cel_file=ARG[i+1]
		}else if (ARG[i] == "-i") {
			cel_file2=ARG[i+1]
		} else if (ARG[i] == "-o") {
			out_file=ARG[i+1]
		} else if (ARG[i] == "-h") {
			usageSampleQC("")
		}
	}
	if (!(file.exists(cel_file))) {
		usagePlateQC("Error : CEL file list file not found")
	}
	if (!(file.exists(cel_file2))) {
		usagePlateQC("Error : CEL file list for DQC filtering not found")
	}
	if (!(file.exists(dqc_file))) {
		usageSampleQC("Error : DQC rate file not found")
	}
	if (!(file.exists(cr_file))) {
		usageSampleQC("Error : genotype call rate file not found")
	}
	if (!(file.exists(match_file))) {
		usagePlateQC("Error : cell to plate match file not found")
	}
	if (!(is.numeric(min_CR))) {
		usagePlateQC("Error : threshold min plate QC should be numeric [0-100]")
	}
	if (!(is.numeric(min_PPR))) {
		usagePlateQC("Error : threshold min pass plate rate should be numeric [0-100]")
	}
	if (out_file == "" ) {
		usagePlateQC("Error : output CEL file list not specified")
	}
	library(pheatmap)
	library(ggplot2)
	library(RColorBrewer)
	library(Cairo)
	celList=read.table(cel_file,header=T)
	cel_table=data.frame(Plate=dirname(as.vector(celList$cel_files)),Sample=basename(as.vector(celList$cel_files)) )
	celList2=read.table(cel_file2,header=T)
	cel_table2=data.frame(Plate=dirname(as.vector(celList2$cel_files)),Sample=basename(as.vector(celList2$cel_files)) )
	match_table=read.table(match_file,header=T)
	cel_table=match_table[match_table$Sample %in% basename(as.vector(celList$cel_files)),]
	cel_table$Plate=as.factor(as.vector(cel_table$Plate))
	cr_table=read.table(cr_file,header=T)
	dqc_table_tmp=read.table(dqc_file,header=T)
	dqc_table=dqc_table_tmp[dqc_table_tmp$cel_files %in% cr_table$cel_files,]

	col=rep(0,dim(cr_table)[1])
	shap=rep(0,dim(cr_table)[1])
	ct=1
	for (i in levels(cr_table$computed_gender)) {
	  shap[cr_table$computed_gender == i]=ct
	  ct=ct+1
	}
	plate_to_keep=NULL
	sampleQC_metrics=data.frame(Plate_Barcode=levels(cel_table$Plate),Result=".",Initial_Sample_Number=0,sample_failing_DQC=0,sample_failing_QC_CR=0,sample_pass=0,percent_passed=0,average_CR_passed_sample=0)
	#jpeg(paste(dirname(cr_file),"Qc_Call_Rate_byPlate.jpg",sep="/"))
	#layout(matrix(1:length(levels(cel_table$Plate)),ncol=1))
	ct=1
	for (i in levels(cel_table$Plate)) {
		print(paste("ploting plate",i))
		CairoJPEG(filename=paste(dirname(cr_file),paste("Qc_Call_Rate_byPlate",i,"jpg",sep="."),sep="/"),width=800,height=400)
		col[dqc_table$cel_files %in% as.vector(match_table[match_table$Plate == i,2])]=ct
		df = data.frame(call_rate=cr_table$call_rate[cr_table$cel_files %in% as.vector(match_table[match_table$Plate == i,2])],affymetrix.plate.peg.wellposition=cr_table$affymetrix.plate.peg.wellposition[cr_table$cel_files %in% as.vector(match_table[match_table$Plate == i,2])])
		## Split column
		df$y = gsub("(.)..", "\\1", df$affymetrix.plate.peg.wellposition)
		df$x = gsub(".(..)", "\\1", df$affymetrix.plate.peg.wellposition)
		## Reverse Y label order
		df$y = factor(df$y, levels=sort(unique(df$y), decreasing=TRUE))
		## Plot
		pal = brewer.pal(n = 10, name =  "RdYlBu")
		#pal = c(rep(pal[1],4), pal)
		p = ggplot(df, aes(x=x, y=y, fill=call_rate)) + geom_tile() + scale_fill_gradientn(name="Call Rate", colours = pal,limits=c(50,100),values=c(0,0.92,seq(0.93,1,length.out=8))) + geom_text(aes(label=round(call_rate,2))) + ggtitle(i) + xlab("") + ylab("")
		print(p)
		dev.off()
		total_cel_num=dim(match_table[match_table$Plate == i ,])[1]
		sampleQC_metrics$Initial_Sample_Number[sampleQC_metrics$Plate_Barcode == i]=total_cel_num
		filtered_cel_num=dim(cel_table[cel_table$Plate == i ,])[1]
		filtered_cel_num2=dim(cel_table2[cel_table2$Plate == i ,])[1]
		sampleQC_metrics$sample_failing_DQC[sampleQC_metrics$Plate_Barcode == i]=total_cel_num-filtered_cel_num2
		sampleQC_metrics$sample_failing_QC_CR[sampleQC_metrics$Plate_Barcode == i]=filtered_cel_num-filtered_cel_num2
		sampleQC_metrics$sample_pass[sampleQC_metrics$Plate_Barcode == i]=filtered_cel_num
		plate_pass_rate=100*(filtered_cel_num/total_cel_num)
		sampleQC_metrics$percent_passed[sampleQC_metrics$Plate_Barcode == i]=plate_pass_rate
		Average_call_rate=mean(cr_table$call_rate[cr_table$cel_files %in% as.vector(match_table[match_table$Plate == i ,2])])
		sampleQC_metrics$average_CR_passed_sample[sampleQC_metrics$Plate_Barcode == i]=
		if (plate_pass_rate > min_PPR) {
			if (Average_call_rate > min_CR) {
				plate_to_keep=c(plate_to_keep,i)
				sampleQC_metrics$Result[sampleQC_metrics$Plate_Barcode == i]="FAILED"
			} else {
				print(paste("Warning low average call rate (",as.character(Average_call_rate),") - Excluding plate",i,"\n",sep=" "))
				sampleQC_metrics$Result[sampleQC_metrics$Plate_Barcode == i]="PASSED"
			}
		} else {
			print(paste("Warning low pass sample rate (",as.character(plate_pass_rate),")  - Excluding plate",i,"\n",sep=" "))
			sampleQC_metrics$Result[sampleQC_metrics$Plate_Barcode == i]="PASSED"
		}
		sample_to_keep=cr_table$cel_files[cr_table$call_rate >= min_CR]
		ct=ct+1
        }
        jpeg(paste(dirname(cr_file),"Qc_Call_Rate_vs_DQC.jpg",sep="/"),res=300)
        par(xpd=TRUE)
        par(mar=c(5, 4, 4, 12))
        plot(x=dqc_table$axiom_dishqc_DQC[match(dqc_table$cel_files,cr_table$cel_files)],y=cr_table$call_rate,col=col,pch=15+shap[match(dqc_table$cel_files,cr_table$cel_files)],xlab="Dish QC",ylab="QC Call Rate",main="QC Call Rate vs. Dish QC")
        legend(1.01,100,legend=levels(cel_table$Plate),fill=unique(sort(col)))
        legend(1.01,98,legend=levels(cr_table$computed_gender),pch=16:(15+length(levels(cr_table$computed_gender))))
        dev.off()
	celList_filtered=data.frame(cel_files=celList$cel_files[dirname(as.vector(celList$cel_files)) %in% plate_to_keep])
	write.table(celList_filtered,out_file,col.names=T,row.names=F,quote=F)
	write.table(sampleQC_metrics,paste(dirname(cr_file),"Plate_QC_metrics.tsv",sep="/"),col.names=T,row.names=F,quote=F,sep="\t")
}

plotSNPmain=function(ARG) {
	#get specific args
	snp_number=6
	output_dir=""
	for (i in 1:length(ARG)) {
		if (ARG[i] == "-d") {
			snp_number=as.integer(ARG[i+1])
		} else if (ARG[i] == "-o") {
			output_dir=ARG[i+1]
		} else if (ARG[i] == "-h") {
			usagePlotSNP("")
		}
	}
	if (!(is.integer(snp_number))) {
		usagePlotSNP("Error : mumber of SNP to plot must be an integer")
	}
	if (output_dir == "" ) {
		usagePlotSNP("Error : output folder not specified")
	}
	library("SNPolisher")
	snp_type=c("CallRateBelowThreshold","Hemizygous","MonoHighResolution","NoMinorHom","OffTargetVariant","Other","PolyHighResolution","Recommended")
	dir.create("tempDir", showWarnings=F)

	for (i in snp_type) {
		Ps_Visualization(pidFile=paste(paste(output_dir,i,sep="/"),"ps",sep="."), output.pdfFile=paste(output_dir,"/Cluster_",i,".pdf",sep=""), summaryFile=paste(output_dir,"AxiomGT1.summary.txt",sep="/"), callFile=paste(output_dir,"AxiomGT1.calls.txt",sep="/"), confidenceFile=paste(output_dir,"AxiomGT1.confidences.txt",sep="/"), posteriorFile=paste(output_dir,"AxiomGT1.snp-posteriors.txt",sep="/"), temp.dir=paste(output_dir,"tempDir/",sep="/"), refFile=NULL, plot.prior=T, priorFile=NULL, match.cel.file.name=F, max.num.SNP.draw=snp_number)
	}

}

mergeOutput=function(ARG) {
	#get specific args
	chunk_number=as.integer(10)
	output_dir=""
	for (i in 1:length(ARG)) {
		if (ARG[i] == "-n") {
			chunk_number=as.integer(ARG[i+1])
		} else if (ARG[i] == "-o") {
			output_dir=ARG[i+1]
		} else if (ARG[i] == "-h") {
			usageMerge("")
		}
	}
	if (!(is.integer(chunk_number))) {
		usageMerge("Error : mumber of chunks must be an integer")
	}
	if (output_dir == "" ) {
		usageMerge("Error : output folder not specified")
	}
	fi=c("AxiomGT1.snp-posteriors.txt","AxiomGT1.calls.txt","AxiomGT1.summary.txt","AxiomGT1.confidences.txt")
	for (i in fi) {
	  if (file.exists(paste(output_dir,"tmp1",i,sep="/"))){
	    print(paste("Merging",i,"chunks"))
	    tmpH=readLines(paste(output_dir,"tmp1",i,sep="/"))
	    headerC=tmpH[grep("^#",tmpH)]
	    write(headerC,file=paste(output_dir,i,sep="/"))
	    for  (j in 1:10) {
	      data=read.table(paste(output_dir,paste("tmp",as.character(j),sep=""),i,sep="/"),header=T,sep="\t")
	      if (j == 1) {
		write.table(data,file=paste(output_dir,i,sep="/"),quote=F,col.names=T,row.names=F,sep="\t",append=T)
	      } else {
		write.table(data,file=paste(output_dir,i,sep="/"),quote=F,col.names=F,row.names=F,sep="\t",append=T)
	      }
	    }
	  }
	}
	fi=c("AxiomGT1.report.txt")
	for (i in fi) {
	  if (file.exists(paste(output_dir,"tmp1",i,sep="/"))){
	    print(paste("combining",i,"chunks"))
	    tmpH=readLines(paste(output_dir,"tmp1",i,sep="/"))
	    headerC=tmpH[grep("^#",tmpH)]
	    write(headerC,file=paste(output_dir,i,sep="/"))
	    for  (j in 1:10) {
	      if (j == 1) {
		data=read.table(paste(output_dir,paste("tmp",as.character(j),sep=""),i,sep="/"),header=T,sep="\t")
	      } else {
		datatmp=read.table(paste(output_dir,paste("tmp",as.character(j),sep=""),i,sep="/"),header=T,sep="\t")
		for (k in 1:dim(data)[2]) {
		  if (is.numeric(data[k])) {
		    data[k]=data[k]+datatmp[k]
		  }
		}
	      }
	    }
	    for (k in 1:dim(data)[2]) {
	      if (is.numeric(data[k])) {
		data[k]=data[k]/10
	      }
	    }
	  }
	  write.table(data,file=paste(output_dir,i,sep="/"),quote=F,col.names=T,row.names=F,sep="\t",append=T)
	}
}

# ## check arg consitency
# if (!(file.exists(design_file))) {
# 	usage("Error : Design file not found")
# }

main()

