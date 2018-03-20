##!/bin/env sh
# 
##2018-03-15 - Edouard Henrion - edouard.henrion@computationalgenomics.ca
#Usage createBaitMapFile.sh INPUT_FILE($1) TMP_FILE($2) OUTPUT_FILE($3)

INPUT_FILE=$1
TMP_FILE=$2
OUTPUT_FILE=$3

awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {if ($8 == ".") {id = $5":"$6"-"$7} else {id = $8} print $5,$6,$7,id}' $INPUT_FILE > $TMP_FILE && \
awk '!a[$0]++' $TMP_FILE > $OUTPUT_FILE && \
rm $TMP_FILE