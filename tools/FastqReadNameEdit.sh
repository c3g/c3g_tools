#!/bin/env bash

# Removes the added /1 and /2 by picard's sam_to_fastq transformation to avoid issues with downstream software like HOMER
# Originaly mplemented by Rola Dali within the hicseq.py script
# Re-implemented as a MUGQIC_TOOLS by Edouard Henrion - 28/02/2018

usage() { 
  echo "Usage: FastqReadNameEdit.sh -i <path> -o <path> -p <path>"
  echo "          [-i <Path of the input fastq.gz file>]"
  echo "          [-o <Path of the output fastq.gz file>]" 
  echo "          [-p <Absolute path of the input fastq.gz>]" 

  exit 1
}

INPUT_FASTQ=""
OUTPUT_FASTQ=""
FASTQ_ABS_PATH=""

while getopts ":i:o:p:" o; do
    case "${o}" in
        i)
            INPUT_FASTQ=${OPTARG}
            ;;
        o)
            OUTPUT_FASTQ=${OPTARG}
            ;;
        p)
            FASTQ_ABS_PATH=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

if [[ ! -s $INPUT_FASTQ ]]
then
  echo "ERROR: The ${INPUT_FASTQ} file doesn't exist or is empty." >&2
  exit 1;
fi

if [[ ! -s $FASTQ_ABS_PATH ]]
then
  echo "ERROR: The ${FASTQ_ABS_PATH} file doesn't exist or is empty." >&2
  exit 1;
fi


## assumes reads in fastq file start with @; if not change
readID=$(zcat $INPUT_FASTQ | head -n 1)
if grep -q '^@.*/[12]$' <<< $readID; then
  zcat $INPUT_FASTQ | sed '/^@/s/\/[12]\>//g' | gzip > $OUTPUT_FASTQ
else
  ln -s -f $FASTQ_ABS_PATH $OUTPUT_FASTQ
fi