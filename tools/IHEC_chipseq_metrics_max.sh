#!/bin/bash


# Computation of IHEC ChIP-seq quality metrics
# Originally from https://github.com/IHEC/ihec-assay-standards/tree/master/ChIP-seq_QC
# reimplemented by mathieu.bourgey@mcgill.ca - 31/07/2017 
# module to load samtools mugqic_dev/deeptools/2.5.3 

usage() { 
  echo "Usage: IHEC_chipseq_metrics.sh [option] [-t narrow|broad|Input]"
  echo "          [-d <ChIP markDup bam>]"
  echo "          [-i <Input markDup bam]" 
  echo "          [-s <ChIp Sample name]" 
  echo "          [-j <Input sample name]" 
  echo "          [-o <Output directory]" 
  echo "          [-p <ChIP_bed_file>]" 
  echo "          [-n <threads>]"
  echo "          [-a <assembly>]"
  exit 1
 }

## By default the number of threads is set to 1;
n=1
CHIP_TYPE=""
CHIP_BAM=""
INPUT_BAM=""
INPUT_NAME=""
CHIP_BED_FILE=""
SAMPLE_NAME=""
OUTPUT_DIR="ihec_metrics"

while getopts "t:d:i:j:p:s:o:n:a::" o; do
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
        j)
            INPUT_NAME=${OPTARG}
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
        a)
            assembly=${OPTARG}
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

if [[ ! ( "${CHIP_TYPE}" == "narrow" || "${CHIP_TYPE}" == "broad" || "${CHIP_TYPE}" == "Input" ) ]]
then
  echo "The experiment type defined isn't one of the following." >&2
  echo "narrow | broad | Input" >&2
  exit 1;
fi 


## need to run script with samples without an input
# if [[  ! -s $INPUT_BAM ]]
# then
#   echo "ERROR: File ${INPUT_BAM} doesn't exist or is empty. The matched Input sample you provided doesn't seem to have been preprocessed yet." >&2
#   exit 1;
# fi

if [ $INPUT_NAME == "no_input" ]
then
  echo "... sample has no input ..."
	INPUT_BAM=""
fi	


## The original number of reads and the number of those aligned:
samtools flagstat ${CHIP_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt

raw_reads_chip=$(awk -v SAMPLE_NAME=$SAMPLE_NAME '{if ($1 == SAMPLE_NAME) print $0}' metrics/trimSampleTable.tsv | cut -f 2)
total_reads_chip=`grep "in total" ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* in total .*//'`
mapped_reads_chip=`grep "mapped (" ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
dupped_reads_chip=`grep "duplicates" ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* duplicates$//'`
dup_rate_chip=$(echo "${dupped_reads_chip}/${mapped_reads_chip}" | bc -l)
aln_rate_chip=$(echo "${mapped_reads_chip}/${total_reads_chip}" | bc -l)
filt_rate_chip=$(echo "${total_reads_chip}/${raw_reads_chip}" | bc -l)

## Finally, the number of singletons for paired-end data sets can be calculated using:
singletons_chip=`grep "singletons" ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* singletons .*//'`


## Remove unmapped read, duplicate reads and those with mapping quality less than 5:
samtools view -b -F 3844 -q 5  ${CHIP_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam

## Index the final deduplicated BAM file
samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam

## run on the input if provided:
if [[ -s $INPUT_BAM ]]
then
  samtools flagstat ${INPUT_BAM} > ${OUTPUT_DIR}/${INPUT_NAME}.markDup_flagstat.txt

  raw_reads_input=$(awk -v INPUT_NAME=$INPUT_NAME '{if ($1 == INPUT_NAME) print $0}' metrics/trimSampleTable.tsv | cut -f 2)
  total_reads_input=`grep "in total" ${OUTPUT_DIR}/${INPUT_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* in total .*//'`
  mapped_reads_input=`grep "mapped (" ${OUTPUT_DIR}/${INPUT_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
  dupped_reads_input=`grep "duplicates" ${OUTPUT_DIR}/${INPUT_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* duplicates$//'`
  dup_rate_input=$(echo "${dupped_reads_input}/${mapped_reads_input}" | bc -l)
  aln_rate_input=$(echo "${mapped_reads_input}/${total_reads_input}" | bc -l)
  filt_rate_input=$(echo "${total_reads_input}/${raw_reads_input}" | bc -l)

  ## Finally, the number of singletons for paired-end data sets can be calculated using:
  singletons_input=`grep "singletons" ${OUTPUT_DIR}/${INPUT_NAME}.markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* singletons .*//'`

  samtools view -b -F 3844 -q 5  ${INPUT_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam
  samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam
  final_reads_input=`samtools flagstat  ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam | grep "mapped (" | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
fi


final_reads_chip=`samtools flagstat  ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam | grep "mapped (" | sed -e 's/ + [[:digit:]]* mapped (.*)//'`


#3.     Calculating Jensen-Shannon distance (JSD)

#To calculate the Jensen-Shannon distance we run:
## Attention: Regarding the bin size (specified in the command below by the ‘-bs’ option) there hasn’t been an agreement on what the optimal bin size is yet. There have been discussions on adopting smaller bin sizes for the sharp peaks and larger bin sizes for the broad peaks.
## No need to remove the blacklisted regions for the JSD calculation.

if [[ "${CHIP_TYPE}" == "narrow" ]]
then
  bin_size=200
else
  bin_size=1000
fi

echo "Experiment type: ${CHIP_TYPE} and bin size: $bin_size" >&2

if [[ -s $INPUT_BAM ]]
then

plotFingerprint -b ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam -bs ${bin_size} -l ${SAMPLE_NAME} INPUT_${CHIP_TYPE} --JSDsample ${OUTPUT_DIR}/${SAMPLE_NAME}_IMPUT.dedup.bam --outQualityMetrics ${OUTPUT_DIR}/${SAMPLE_NAME}.fingerprint.txt -plot ${OUTPUT_DIR}/${SAMPLE_NAME}.fingerprint.png -p $n

js_dist=`grep ${SAMPLE_NAME} ${OUTPUT_DIR}/${SAMPLE_NAME}.fingerprint.txt | cut -f 8`
chance_div=`grep ${SAMPLE_NAME} ${OUTPUT_DIR}/${SAMPLE_NAME}.fingerprint.txt | cut -f 12`

fi

#4.     Calculating FRiP scores
nmb_peaks=$(wc -l ${CHIP_BED_FILE} | cut -f 1 -d " ")
reads_under_peaks=`samtools view -c -L ${CHIP_BED_FILE} ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam`
frip=$(echo "${reads_under_peaks}/${final_reads_chip}" | bc -l)

#5. extract NSC and RSC from run_spp

nsc_chip=$(cut -f 9 ${OUTPUT_DIR}/${SAMPLE_NAME}.crosscor | head -n 1)
rsc_chip=$(cut -f 10 ${OUTPUT_DIR}/${SAMPLE_NAME}.crosscor | head -n 1)
quality_chip_num=$(cut -f 11 ${OUTPUT_DIR}/${SAMPLE_NAME}.crosscor | head -n 1)

## Quality tag based on thresholded RSC (codes= -2:veryLow, -1:Low, 0:Medium, 1:High, 2:veryHigh)

if [[ "$quality_chip_num" == "-2" ]]
  then
    quality_chip=veryLow
elif [[ "$quality_chip_num" == "-1" ]]
  then
    quality_chip=Low
elif [[ "$quality_chip_num" == "0" ]]
  then
    quality_chip=Medium
elif [[ "$quality_chip_num" == "1" ]]
  then
    quality_chip=High
elif [[ "$quality_chip_num" == "2" ]]
  then
    quality_chip=veryHigh
fi





if [[ -s $INPUT_BAM ]]
then
nsc_input=$(cut -f 9 ${OUTPUT_DIR}/${INPUT_NAME}.crosscor | head -n 1)
rsc_input=$(cut -f 10 ${OUTPUT_DIR}/${INPUT_NAME}.crosscor | head -n 1)
quality_input_num=$(cut -f 11 ${OUTPUT_DIR}/${INPUT_NAME}.crosscor | head -n 1)


if [[ "$quality_input_num" == "-2" ]]
  then
    quality_input=veryLow
elif [[ "$quality_input_num" == "-1" ]]
  then
    quality_input=Low
elif [[ "$quality_input_num" == "0" ]]
  then
    quality_input=Medium
elif [[ "$quality_input_num" == "1" ]]
  then
    quality_input=High
elif [[ "$quality_input_num" == "2" ]]
  then
    quality_input=veryHigh
fi



fi


LC_NUMERIC="en_US.UTF-8"

#printf "genome_assembly\ttreat_name\tctl_name\ttreat_filtered_reads\ttreat_mapped_reads\tclt_filtered_reads\tclt_mapped_reads\ttreat_aln_frac\tctl_aln_frac\ttreat_dup_frac\tctl_dup_frac\tnmb_peaks\treads_in_peaks\tfrip\ttreat_nsc\tctrl_nsc\ttreat_rsc\tctrl_rsc\ttreat_Quality\tctrl_Quality\tsingletons\tjs_dist\tchance_div\n" > ${OUTPUT_DIR}/${SAMPLE_NAME}.read_stats.txt
#LANG=C printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t\t%.4f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\t%d\t%.4f\t%.4f\n" "${assembly}" "${SAMPLE_NAME}" "${INPUT_NAME}" "$total_reads_chip" "$mapped_reads_chip" "$total_reads_input" "$mapped_reads_input" "$aln_rate_chip" "$aln_rate_input" "$dup_rate_chip" "$dup_rate_input" "$nmb_peaks" "$reads_under_peaks" "$frip" "$nsc_chip" "$nsc_input" "$rsc_chip" "$rsc_input" "$quality_chip" "$quality_input" "$singletons_chip" "$js_dist" "$chance_div"  >> ${OUTPUT_DIR}/${SAMPLE_NAME}.read_stats.txt

printf "genome_assembly\tChIP_type\ttreat_name\tctl_name\ttreat_raw_reads\ttreat_filtered_reads\ttreat_mapped_reads\ttreat_duplicated_reads\ttreat_final_reads\tctl_raw_reads\tclt_filtered_reads\tclt_mapped_reads\tctl_duplicated_reads\tctl_final_reads\ttreat_filtered_frac\tctl_filtered_frac\ttreat_aln_frac\tctl_aln_frac\ttreat_dup_frac\tctl_dup_frac\tnmb_peaks\treads_in_peaks\tfrip\ttreat_nsc\tctrl_nsc\ttreat_rsc\tctrl_rsc\ttreat_Quality\tctrl_Quality\tsingletons\tjs_dist\tchance_div\n" > ${OUTPUT_DIR}/IHEC_metrics_chipseq_${SAMPLE_NAME}.txt
LANG=C printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\t%d\t%.4f\t%.4f\n" "${assembly}" "${CHIP_TYPE}" "${SAMPLE_NAME}" "${INPUT_NAME}" "$raw_reads_chip" "$total_reads_chip" "$mapped_reads_chip" "$dupped_reads_chip" "$final_reads_chip" "$raw_reads_input" "$total_reads_input" "$mapped_reads_input" "$dupped_reads_input" "$final_reads_input" "$filt_rate_chip" "$filt_rate_input" "$aln_rate_chip" "$aln_rate_input" "$dup_rate_chip" "$dup_rate_input" "$nmb_peaks" "$reads_under_peaks" "$frip" "$nsc_chip" "$nsc_input" "$rsc_chip" "$rsc_input" "$quality_chip" "$quality_input" "$singletons_chip" "$js_dist" "$chance_div"  >> ${OUTPUT_DIR}/IHEC_metrics_chipseq_${SAMPLE_NAME}.txt




