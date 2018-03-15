#!/bin/env bash

if [ $# -ne 3 ]; then
    echo $0: usage: CreateHicFile.sh myNameSortedbam sampleId tmpDir
    exit 1
fi
bam=$1
id=$2
tmpDir=$3

## check if bam is sorted by name, otherwise abort:
## check if bam is sorted by name, otherwise abort:
sort_flag=$(samtools view -H ${bam} | grep -oP 'SO:.*')

if [ "$sort_flag" != "SO:queryname" ]; then
    echo "input bam must be sorted by query name! please sort bam!"
    exit 1
fi

##extract read info to make .hic file
## remove chrs containing random, chrUn, hap, chrM, chrY

echo "... ... ... creating ${id}.juicebox.input"

samtools view ${bam} | awk 'BEGIN {FS="\t"; OFS="\t"} {name1=$1; str1=and($2,16); chr1=substr($3, 4); pos1=$4; mapq1=$5; getline; name2=$1; str2=and($2,16); chr2=substr($3, 4); pos2=$4; mapq2=$5; 
	if (chr1 == "X"){chr1=23}
	if (chr2 == "X"){chr2=23}
	if (name1==name2) {
		#print chr1, chr2
		if ((chr1 ~ /random|Un|hap|M|Y/) || (chr2 ~ /random|Un|hap|M|Y/)) {
		#print "match found!"
		}
		else {
		#print "entered inner loop"
			if (chr1>chr2) {print NR, str2, chr2, pos2,1, str1, chr1, pos1, 0, mapq2, mapq1} 
			else {print NR, str1, chr1, pos1, 0, str2, chr2, pos2 ,1, mapq1, mapq2}	
		}
	}
	else {getline}}' > ${id}.juicebox.input

## sort file by columns 3 and 7:
## needs big scratch directory

echo "... ... ... sorting ${id}.juicebox.input"

sort -T $tmpDir -k3,3d -k7,7d < ${id}.juicebox.input > ${id}.juicebox.input.sorted


