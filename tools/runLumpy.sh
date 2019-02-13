#!/bin/env sh
set -e
set -o pipefail

export LUMPY_SV_HOME=/lb/project/mugqic/analyste_dev/software/lumpy-sv/lumpy-sv-0.2.9/
export PATH=${LUMPY_SV_HOME}/bin:$PATH


BLOOD=$1
TUMOR=$2
BLOOD_PE=$3
TUMOR_PE=$4
BLOOD_SR=$5
TUMOR_SR=$6
READ_LENGTH=$7
OUTPUT=$8

echo "Extracting iSize histograms " `date +"%Y-%m-%d %T"`
# Can't use pipefail because some pipes get closed too early by design
set +o pipefail
samtools view $BLOOD | tail -n+100000 | python $LUMPY_SV_HOME/scripts/pairend_distro.py -r $READ_LENGTH -X 4 -N 10000 -o ${BLOOD_PE%.bam}.histo > ${BLOOD_PE%.bam}.histo.out
samtools view $TUMOR | tail -n+100000 | python $LUMPY_SV_HOME/scripts/pairend_distro.py -r $READ_LENGTH -X 4 -N 10000 -o ${TUMOR_PE%.bam}.histo > ${TUMOR_PE%.bam}.histo.out
bloodDist=$(cat ${BLOOD_PE%.bam}.histo.out)
tumorDist=$(cat ${TUMOR_PE%.bam}.histo.out)

bloodMean=`echo $bloodDist | sed 's/mean:\([^ ]\+\) .*/\1/g'`
bloodStdev=`echo $bloodDist | sed 's/.* stdev:\([0-9.]\+\).*/\1/g'`

tumorMean=`echo $tumorDist | sed 's/mean:\([^ ]\+\) .*/\1/g'`
tumorStdev=`echo $tumorDist | sed 's/.* stdev:\([0-9.]\+\).*/\1/g'`

echo "Running  lumpy-sv" `date +"%Y-%m-%d %T"`
lumpy -mw 4 -tt 0.0 \
  -pe bam_file:${TUMOR_PE},histo_file:${TUMOR_PE%.bam}.histo,mean:$tumorMean,stdev:$tumorStdev,read_length:$READ_LENGTH,min_non_overlap:$READ_LENGTH,discordant_z:4,back_distance:20,weight:1,id:10,min_mapping_threshold:1\
  -pe bam_file:${BLOOD_PE},histo_file:${BLOOD_PE%.bam}.histo,mean:$bloodMean,stdev:$bloodStdev,read_length:$READ_LENGTH,min_non_overlap:$READ_LENGTH,discordant_z:4,back_distance:20,weight:1,id:20,min_mapping_threshold:1 \
  -sr bam_file:${TUMOR_SR},back_distance:20,weight:1,id:11,min_mapping_threshold:20 \
  -sr bam_file:${BLOOD_SR},back_distance:20,weight:1,id:21,min_mapping_threshold:20 \
  > $OUTPUT

echo "Done " `date +"%Y-%m-%d %T"`
