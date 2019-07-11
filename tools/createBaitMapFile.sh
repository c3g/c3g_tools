#!/usr/bin/env sh
# 
##2018-03-15 - Edouard Henrion - edouard.henrion@computationalgenomics.ca
#Usage createBaitMapFile.sh INPUT_FILE($1) BAIT_FILE($2) SORTED_BAIT_FILE($3) ANNOTATION($4) TMP_FILE($5) OUTPUT_FILE($6)

INPUT_FILE=$1
BAIT_FILE=$2
SORTED_BAIT_FILE=$3
ANNOTATION=$4
TMP_FILE=$5
OUTPUT_FILE=$6

column_num=$(awk 'NR <2 {{print NF}}' $BAIT_FILE)

## annotate file with annotation in baitBed otherwise annotate with random id
if [[ $column_num -eq 4 ]]; then
    bedmap --echo --echo-map-id $INPUT_FILE $SORTED_BAIT_FILE \
    | tr '|' '\t' > $OUTPUT_FILE;
else
    awk -v ANNOTATION="$ANNOTATION" 'BEGIN {{FS="\t"; OFS="\t"}}{{print $0, "ANNOTATION"NR}}' $INPUT_FILE > $OUTPUT_FILE
fi
