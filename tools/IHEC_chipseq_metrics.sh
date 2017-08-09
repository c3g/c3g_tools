#!/bin/bash


# Computation of IHEC ChIP-seq quality metrics
# Originally from https://github.com/IHEC/ihec-assay-standards/tree/master/ChIP-seq_QC
# reimplemented by mathieu.bourgey@mcgill.ca - 31/07/2017 
# mv samtools to sambamba
# module to load sambamba picard mugqic_dev/deeptools/2.5.3 (for deeptools)
usage() { 
  echo "Usage: IHEC_chipseq_metrics.sh [option] [-t <H3K27ac|H3K27me3|H3K36me3|H3K4me1|H3K4me3|H3K9me3|Input|H2AFZ|H3ac|H3K4me2|H3K9ac>]"
  echo "          [-d <ChIP_file_prefix>]"
  echo "          [-u <Input_file_prefix]" 
  echo "          [-p <ChIP_bed_file>]" 
  echo "          [-n <threads>]"
  exit 1
 }

## By default the number of threads is set to 1;
n=1
DEDUP_BAM=""




while getopts "t:d:i:p:s:o:n::" o; do
    case "${o}" in
        t)
            CHIP_TYPE=${OPTARG}
            ;;
        d)
            CHIP_BAM=${OPTARG}
            ;;
        i)
            INPUT_BAM=${OPTARG}
            ;;
        p)
            CHIP_BED_FILE=${OPTARG}
            ;;
        s)
            SAMPLE_NAME=${OPTARG}
            ;;
        o)
            OUTPUT_DIR=${OPTARG}
            ;;
        n)
            n=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done


## If the ChIP bam file doesn't exist, then throw an error 
if [[  ! -s $CHIP_BAM ]]
then
  echo "ERROR: The ${CHIP_BAM} file doesn't exist or is empty." >&2
  exit 1;
fi

if [[ ( -z "${CHIP_BED_FILE}" || ! -s ${CHIP_BED_FILE} ) && ! "${CHIP_TYPE}" == "Input" ]]
then
  echo "ERROR: Your sample isn't of type Input but you're not providing a BED file with the peaks, or the file doesn't exist or is empty." >&2
  exit 1;
fi

if [[ -z "${CHIP_TYPE}" || -z "${SAMPLE_NAME}"  ]]
then
  usage
fi

if [[ ! ( "${CHIP_TYPE}" == "H3K27ac" || "${CHIP_TYPE}" == "H3K27me3" || "${CHIP_TYPE}" == "H3K36me3" || "${CHIP_TYPE}" == "H3K4me1" || "${CHIP_TYPE}" == "H3K4me3" || "${CHIP_TYPE}" == "H3K9me3" || "${CHIP_TYPE}" == "Input" || "${CHIP_TYPE}" == "H2AFZ" || "${CHIP_TYPE}" == "H3ac" || "${CHIP_TYPE}" == "H3K4me2" || "${CHIP_TYPE}" == "H3K9ac" ) ]]
then
  echo "The experiment type defined isn't one of the following." >&2
  echo "H3K27ac | H3K27me3 | H3K36me3 | H3K4me1 | H3K4me3 | H3K9me3 | Input | H2AFZ | H3ac | H3K4me2 | H3K9ac" >&2
  exit 1;
fi 


if [[  ! -s $INPUT_BAM ]]
then
  echo "ERROR: File ${INPUT_BAM} doesn't exist or is empty. The matched Input sample you provided doesn't seem to have been preprocessed yet." >&2
  exit 1;
fi


## The original number of reads and the number of those aligned:
sambamba flagstat ${CHIP_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt


total_reads=`grep "in total" ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* in total .*//'`
mapped_reads=`grep "mapped (" ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
dupped_reads=`grep "duplicates" ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* duplicates$//'`
dup_rate=$(echo "${dupped_reads}/${mapped_reads}" | bc -l)

## Finally, the number of singletons for paired-end data sets can be calculated using:
singletons=`grep "singletons" ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* singletons .*//'`


## Remove unmapped read, duplicate reads and those with mapping quality less than 5:
sambamba view -b -F 3844 -q 5  ${CHIP_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam
sambamba view -b -F 3844 -q 5  ${INPUT_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam

## Index the final deduplicated BAM file
sambamba index ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam
sambamba index ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam


final_reads=`sambamba flagstat  ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam | grep "mapped (" | sed -e 's/ + [[:digit:]]* mapped (.*)//'`

#3.     Calculating Jensen-Shannon distance (JSD)

#To calculate the Jensen-Shannon distance we run:
## Attention: Regarding the bin size (specified in the command below by the ‘-bs’ option) there hasn’t been an agreement on what the optimal bin size is yet. There have been discussions on adopting smaller bin sizes for the sharp peaks and larger bin sizes for the broad peaks.
## No need to remove the blacklisted regions for the JSD calculation.

if [[ "${CHIP_TYPE}" == "H3K27ac" || "${CHIP_TYPE}" == "H3K4me3" || "${CHIP_TYPE}" == "H2AFZ" || "${CHIP_TYPE}" == "H3ac" || "${CHIP_TYPE}" == "H3K4me2" || "${CHIP_TYPE}" == "H3K9ac" ]]
then
  bin_size=200
else
  bin_size=1000
fi

echo "Experiment type: ${CHIP_TYPE} and bin size: $bin_size" >&2

plotFingerprint -b ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam -bs ${bin_size} -l ${SAMPLE_NAME} INPUT_${CHIP_TYPE} --JSDsample ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam --outQualityMetrics ${OUTPUT_DIR}/${SAMPLE_NAME}.fingerprint.txt -plot ${OUTPUT_DIR}/${SAMPLE_NAME}.fingerprint.png -p $n

js_dist=`grep ${cname} ${OUTPUT_DIR}/${SAMPLE_NAME}.fingerprint.txt | cut -f 8`
chance_div=`grep ${cname} ${OUTPUT_DIR}/${SAMPLE_NAME}.fingerprint.txt | cut -f 12`

#4.     Calculating FRiP scores
reads_under_peaks=`sambamba view -c -L ${CHIP_BED_FILE} ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam`
frip=$(echo "${reads_under_peaks}/${final_reads}" | bc -l)


printf "ChIP_name\tInput_name\ttotal_reads\tmapped_reads\tdupped_reads\tdup_rate\tsingletons\tfinal_reads\tjs_dist\tchance_div\tfrip\n" > $WORKING_DIR/${cname}/${cname}_read_stats.txt
printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\n" "$cname" "$iname" "$total_reads" "$mapped_reads" "$dupped_reads" "$dup_rate" "$singletons" "$final_reads" "$js_dist" "$chance_div" "$frip" >> ${OUTPUT_DIR}/${SAMPLE_NAME}.read_stats.txt

