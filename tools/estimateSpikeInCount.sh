#!/bin/sh
#modified version of runBlast.sh
#by Eloi Mercier

set -e
echo "STARTING: " $(date)

#get args
SAMPLE_TO=$1
OUTPUT_DIR=$2
FILE_1=$3
DATABASE=$4
ALIGN_SIZE=$5
MISTMATCHES=$6
IDENTITY=$7

FILE_NAME=$(basename $FILE_1)

#params
if [[ $FILE_1 =~ \.gz$ ]]; then
	FPR="fastqPickRandom.pl --compressed --threshold "
	read_method=zcat
else
	FPR="fastqPickRandom.pl --threshold "
	read_method=cat
fi
FQ2FA="fastq2FastaQual.pl"
BLAST="blastn -query"

#check output directory
if [ ! -d $OUTPUT_DIR/ ]; then
  mkdir $OUTPUT_DIR
fi

#check database
if [ ! -e $DATABASE.nsd ]; then #.nsd, .nsi, .nsq, .nsh, .nin all created by makeblastdb; check for all?
	echo "generating blast database"
	makeblastdb -in $DATABASE -parse_seqids -dbtype nucl
fi

Nseq=$($read_method $FILE_1 | awk ' { if (substr($0,0,1) == "+") { print $0} }' | wc -l)
echo "$FILE_1 has $Nseq sequences"
#Get the threshold of random picking
thrC=$(echo " scale=6; $SAMPLE_TO / $Nseq" | bc)
echo "$FILE_1 has 0$thrC rdp threshold"

###################
#1. Random pick
###################
$FPR 0$thrC --input1 $FILE_1 --out1 $OUTPUT_DIR/$FILE_NAME.subSampled_${SAMPLE_TO}.fastq

###################
#2. Format fastq to fasta
###################
$FQ2FA $OUTPUT_DIR/$FILE_NAME.subSampled_${SAMPLE_TO}.fastq $OUTPUT_DIR/$FILE_NAME.subSampled_${SAMPLE_TO}.fasta $OUTPUT_DIR/$FILE_NAME.subSampled_${SAMPLE_TO}.qual

###################
#3. Blast fasta
###################
$BLAST $OUTPUT_DIR/$FILE_NAME.subSampled_${SAMPLE_TO}.fasta -db $DATABASE -out $OUTPUT_DIR/$FILE_NAME.subSampled_${SAMPLE_TO}.blast.tsv -max_target_seqs 1  -outfmt "6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore" -ungapped

###################
#4. Filter blast output
###################

python /home/emercier/bitbucket/spikein_dev/mugqic_tools/python-tools/filterBastOutput.py --input_file $OUTPUT_DIR/$FILE_NAME.subSampled_${SAMPLE_TO}.blast.tsv --tag_file $DATABASE --output_dir $OUTPUT_DIR --align_length_min $ALIGN_SIZE --mismatches_max $MISTMATCHES --perc_identity_min $IDENTITY #log file?

###################
#5. Create report and graphs
###################
if [[ $IDENTITY != *.* ]]; then #treat 0 and 1 cases -> 0.0 and 1.1
	IDENTITY=$IDENTITY.0
fi
IDENTITY_PERC=${IDENTITY//./} #0.0 -> 00
########REMOVE HARD CODED PATH
Rscript --vanilla /home/emercier/bitbucket/spikein_dev/mugqic_tools/R-tools/plotSpikeInCount.R $OUTPUT_DIR/$FILE_NAME.subSampled_${SAMPLE_TO}.blast.tsv.${ALIGN_SIZE}bp_${MISTMATCHES}MM_${IDENTITY_PERC}id.tsv $SAMPLE_TO


