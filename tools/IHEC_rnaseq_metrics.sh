#!/bin/sh
set -e
set -o pipefail

# Computation of IHEC RNA-seq quality metrics
# From Sebastain Ullirch <seba.ull@gmail.com>
# mathieu.bourgey@mcgill.ca - 26/07/2017 
# module load sambamba bedtools

SAMPLE_PATH=$1
SAMPLE_NAME=$2
INTERGENIC_BED=$3
RRNA_BED=$4

#create new BAM with only primary mapped and unmapped reads
echo "create new BAM with only primary mapped and unmapped reads..."
sambamba view -b -F 0x900 $SAMPLE_PATH > ${SAMPLE_NAME}_no_multimap.bam
echo "...Done"

#get total number of reads and number of mapped reads
echo "get total number of reads and number of mapped reads"
sambamba flagstat ${SAMPLE_NAME}_no_multimap.bam > ${SAMPLE_NAME}_total_mapped_reads.txt
echo "...Done"

#get the number of reads falling into intergenic regions
echo "get the number of reads falling into intergenic regions"
bedtools coverage -abam ${SAMPLE_NAME}_no_multimap.bam -b $INTERGENIC_BED | awk '{ SUM += $4} END { print SUM }' > ${SAMPLE_NAME}_intergenic_reads.txt
echo "...Done"

#get the number of reads from rRNA
echo "get the number of reads from rRNA"
bedtools coverage -abam ${SAMPLE_NAME}_no_multimap.bam -b $RRNA_BED | awk '{ SUM += $4} END { print SUM }' > ${SAMPLE_NAME}_rRNA_reads.txt
echo "...Done"

#get number of duplicates
echo "get the number of duplicates"
java -Xmx4g -jar picard/MarkDuplicates.jar I=${SAMPLE_NAME}_no_multimap.bam O=${SAMPLE_NAME}_noMULTI_noDUP.bam M=${SAMPLE_NAME}_duplicated.txt
echo "...Done"
