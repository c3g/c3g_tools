# Performs various filter of Axiam gene titan data
# Written by Mathieu Bourgey - August 2016
# Usage : Rscript filterAxiom.R -s path_design -c path_rawcountfile -o output_dir

# Usage functions

usageMain=function(errM) {
	cat("\nUsage : Rscript filterAxiom.R -s <Step>\n")
	cat("       genoqc    : filter sample based on Dish QC value (Best Practices Step 3) \n")
	cat("       sampleqc  : filter sample based on call rate value (Best Practices Step 5) \n")
	cat("       plateqc   : filter sample based on the Plate Pass Rat value (Best Practices Step 6) \n")
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
	cat("       -c        : CEL files list to be filtered\n")
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
	match_file=""
	cr_file=""
	for (i in 1:length(ARG)) {
		if (ARG[i] == "-a") {
			min_CR=as.numeric(ARG[i+1])
		} else if (ARG[i] == "-q") {
			cr_file=ARG[i+1]
		} else if (ARG[i] == "-d") {
			min_PPR=as.numeric(ARG[i+1])
		} else if (ARG[i] == "-p") {
			match_file=ARG[i+1]
		} else if (ARG[i] == "-c") {
			cel_file=ARG[i+1]
		} else if (ARG[i] == "-o") {
			out_file=ARG[i+1]
		} else if (ARG[i] == "-h") {
			usageSampleQC("")
		}
	}
	if (!(file.exists(cel_file))) {
		usagePlateQC("Error : CEL file list file not found")
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
	celList=read.table(cel_file,header=T)
	match_table=read.table(match_file,header=T)
	cel_table=match_table[match_table$Sample %in% basename(as.vector(celList$cel_files)),]
	cr_table=read.table(cr_file,header=T)
	plate_to_keep=NULL
	for (i in levels(cel_table$Plate)) {
		total_cel_num=dim(match_table[match_table$Plate == i ,])[1]
		filtered_cel_num=dim(cel_table[cel_table$Plate == i ,])[1]
		plate_pass_rate=100*(filtered_cel_num/total_cel_num)
		if (plate_pass_rate > min_PPR) {
			Average_call_rate=mean(cr_table$call_rate[cr_table$cel_files %in% as.vector(match_table[match_table$Plate == i ,2])])
			if (Average_call_rate > min_CR) {
				plate_to_keep=c(plate_to_keep,i)
			} else {
				print(paste("Warning low average call rate (",as.character(Average_call_rate),") - Excluding plate",i,"\n",sep=" "))
			}
		} else {
			print(paste("Warning low pass sample rate (",as.character(plate_pass_rate),")  - Excluding plate",i,"\n",sep=" "))
		}
        }
        cf=cel_table[cel_table$Plate %in% plate_to_keep,2]
	celList_filtered=data.frame(cel_files=celList[basename(as.vector(celList$cel_files)) %in% cf,])
	write.table(celList_filtered,out_file,col.names=T,row.names=F,quote=F)
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
# ## check arg consitency
# if (!(file.exists(design_file))) {
# 	usage("Error : Design file not found")
# }

main()

