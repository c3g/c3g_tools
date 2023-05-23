#!/bin/env bash
set -eu -o pipefail

# A Quality Control step for ONT data using BLAST
# Randomly subsamples ONT reads from a dataset
# Trims them to be 1000bp on length
# Then runs BLAST on the trimmed reads and reports the species of the top hit
# Output is a list of all top species hits
# Implemented by H. Galvez using code from other mugqic_tools

if [ "$#" != 3 ] #Are there less/more than three arguments?
then
    echo "Error: you provided an incorrect number of arguments."
    echo "Usage: runBlastQC_ONT.sh output_directory reads_fastq_directory readset_name"
    exit 1
fi


OUTPUT_DIR=$1
READ_DIR=$2
READSET=$3

# Create the output directory (usually called blastqc/SAMPLE)
mkdir -p ${OUTPUT_DIR}

# Using the `file` command, output the file type of the first 5 FASTQ files (used later to check if zipped)
file -L ${READ_DIR}/*.fa* | head >> ${OUTPUT_DIR}/fastq_file_type.tmp

# If the files are zipped, use zcat to generate the input fastq
# If the files are not zipped, use cat to generate input file
if grep -q "gzip" ${OUTPUT_DIR}/fastq_file_type.tmp
 then zcat ${READ_DIR}/*.fastq.gz >> ${OUTPUT_DIR}/full_input.tmp.fastq
elif grep -q "ASCII" ${OUTPUT_DIR}/fastq_file_type.tmp
 then cat ${READ_DIR}/*.fastq >> ${OUTPUT_DIR}/full_input.tmp.fastq
fi

# From input FASTQ file containing all possible reads, count total number of reads
Nseq=$(cat ${OUTPUT_DIR}/full_input.tmp.fastq | awk ' {{ if (substr($0,0,1) == "+") {{ print $0}} }}' | wc -l)

# Then, depending on the total number of reads, decide how many to subsample
thrC=$(echo " scale=6; 1000 / $Nseq" | bc)
if [ $thrC == 0 ]
 then thrC=0.000001
fi

# Use FastqPickRandom.pl from MUGQIC_TOOLS to randomly subsample the determined number of reads
fastqPickRandom.pl --threshold 0$thrC \
  --input1 ${OUTPUT_DIR}/full_input.tmp.fastq \
  --out1 ${OUTPUT_DIR}/subsample_input.fastq

# Remove temporary files
rm ${OUTPUT_DIR}/full_input.tmp.fastq ${OUTPUT_DIR}/fastq_file_type.tmp

# Trim nanopore reads to a length of 1000bp using trim_nanopore.py from MUGQIC_TOOLS
trim_nanopore.py -i ${OUTPUT_DIR}/subsample_input.fastq \
  -o ${OUTPUT_DIR}/subsample_input.trim.fastq -s 1000

# Convert FASTQ to FASTA to prepare for BLAST launch
fastq2FastaQual.pl ${OUTPUT_DIR}/subsample_input.trim.fastq \
  ${OUTPUT_DIR}/subsample_input.trim.fasta \
  ${OUTPUT_DIR}/subsample_input.trim.qual

# Using the FASTA generated above, run a BLASTn query to the nt database
blastn -query ${OUTPUT_DIR}/subsample_input.trim.fasta \
  -db nt -out ${OUTPUT_DIR}/subsample_input.trim.blastres \
  -perc_identity 80 \
  -num_descriptions 1 \
  -num_alignments 1

# From the output of the BLASTn run, parse out the top species hit and output into a summary report file
grep ">" ${OUTPUT_DIR}/subsample_input.trim.blastres | \
  awk ' {{ print $2 "_" $3}} ' | \
  sort | \
  uniq -c | \
  sort -n -r | \
  head -20 > ${OUTPUT_DIR}/${READSET}.blastHit_20MF_species.txt