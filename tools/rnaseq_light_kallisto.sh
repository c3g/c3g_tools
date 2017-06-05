#!/bin/sh
#by Eloi Mercier

set -e
echo "STARTING: " $(date)

#get args
FASTQ1=$1
FASTQ2=$2
TRANSCRIPTOME=$3
GTF=$4
OUTPUT_DIR=$5

# module load mugqic_dev/kallisto/0.43.0
echo $TRANSCRIPTOME
echo $OUTPUT_DIR
echo $FASTQ1
echo $FASTQ2
mkdir -p $OUTPUT_DIR
kallisto quant -i $TRANSCRIPTOME -o $OUTPUT_DIR $FASTQ1 $FASTQ2
mv $OUTPUT_DIR/abundance.tsv $OUTPUT_DIR/abundance_transcripts.tsv

#R script transcript -> gene level
Rscript --vanilla /home/emercier/bitbucket/RNAseq_light_dev/mugqic_tools/R-tools/abundanceTranscript2geneLevel.R $OUTPUT_DIR/abundance_transcripts.tsv $GTF

