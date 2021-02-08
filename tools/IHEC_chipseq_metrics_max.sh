#!/bin/env bash

# Computation of IHEC ChIP-seq quality metrics
# Originally from https://github.com/IHEC/ihec-assay-standards/tree/master/ChIP-seq_QC
# reimplemented by mathieu.bourgey@mcgill.ca - 31/07/2017
# module to load sambamba samtools mugqic/deeptools/2.5.3

usage() { 
  echo "Usage: IHEC_chipseq_metrics.sh [option] [-t narrow|broad|Input]"
  echo "          [-d <ChIP markDup bam>]"
  echo "          [-i <Input markDup bam]"
  echo "          [-s <ChIP Sample name]"
  echo "          [-j <Input Sample name]"
  echo "          [-c <ChIP name]"
  echo "          [-o <Output directory]"
  echo "          [-p <ChIP_bed_file>]"
  echo "          [-n <threads>]"
  echo "          [-a <assembly>]"
  exit 1
 }

## By default the number of threads is set to 1;
n=1
CHIP_NAME=""
CHIP_TYPE=""
CHIP_BAM=""
INPUT_BAM=""
INPUT_NAME=""
CHIP_BED_FILE=""
SAMPLE_NAME=""
OUTPUT_DIR="ihec_metrics"

while getopts "t:d:i:j:c:p:s:o:n:a::" o; do
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
        c)
            CHIP_NAME=${OPTARG}
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
# samtools flagstat ${CHIP_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}.markDup_flagstat.txt
flagstat_file="${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.markDup_flagstat.txt"
sambamba flagstat -t $n ${CHIP_BAM} > $flagstat_file


# raw_reads_chip=$(awk -v SAMPLE_NAME=$SAMPLE_NAME '{if ($1 == SAMPLE_NAME) print $0}' metrics/trimSampleTable.tsv | cut -f 2)
supplementarysecondary_reads_chip=`bc <<< $(grep "secondary" $flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
# trimmed_reads_chip=`bc <<< $(grep "in total" $flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$supplementarysecondary_reads_chip`
mapped_reads_chip=`bc <<< $(grep "mapped (" $flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$supplementarysecondary_reads_chip`
dup_reads_chip=`grep "duplicates" $flagstat_file | sed -e 's/ + [[:digit:]]* duplicates$//'`
dup_rate_chip=$(echo "100*${dup_reads_chip}/${mapped_reads_chip}" | bc -l)
# mapped_rate_chip=$(echo "100*${mapped_reads_chip}/${trimmed_reads_chip}" | bc -l)
# trimmed_rate_chip=$(echo "100*${trimmed_reads_chip}/${raw_reads_chip}" | bc -l)
trimmomatic_table="metrics/trimSampleTable.tsv"
if [[ -s $trimmomatic_table ]]
then
  raw_reads_chip=$(awk -v SAMPLE_NAME=$SAMPLE_NAME '{if ($1 == SAMPLE_NAME) print $0}' $trimmomatic_table | cut -f 2)
  trimmed_reads_chip=`bc <<< $(grep "in total" $flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$supplementarysecondary_reads_chip`
  mapped_rate_chip=$(echo "100*${mapped_reads_chip}/${trimmed_reads_chip}" | bc -l)
  trimmed_rate_chip=$(echo "100*${trimmed_reads_chip}/${raw_reads_chip}" | bc -l)
else
  raw_reads_chip=`bc <<< $(grep "in total" $flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$supplementarysecondary_reads_chip`
  trimmed_reads_chip="NULL"
  mapped_rate_chip=$(echo "100*${mapped_reads_chip}/${raw_reads_chip}" | bc -l)
  trimmed_rate_chip="NULL"
fi


## Finally, the number of singletons for paired-end data sets can be calculated using:
singletons_chip=`grep "singletons" $flagstat_file | sed -e 's/ + [[:digit:]]* singletons .*//'`


## Remove unmapped read, duplicate reads and those with mapping quality less than 5:
# samtools view -b -F 3844 -q 5  ${CHIP_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam
dedup_bam="${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.dedup.bam"
sambamba view -t $n -f bam -F "not unmapped and not secondary_alignment and not failed_quality_control and not duplicate and not supplementary and mapping_quality >= 5" ${CHIP_BAM} > $dedup_bam

## Index the final deduplicated BAM file
# samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam
sambamba index -t $n $dedup_bam

## run on the input if provided:
if [ -s $INPUT_BAM ] && ! [ -s ${OUTPUT_DIR}/${SAMPLE_NAME}.Input.dedup.bam ]
  then
    # samtools flagstat ${INPUT_BAM} > ${OUTPUT_DIR}/${INPUT_NAME}.markDup_flagstat.txt
    input_flagstat_file="${OUTPUT_DIR}/${INPUT_NAME}.Input.markDup_flagstat.txt"
    sambamba flagstat -t $n ${INPUT_BAM} > $input_flagstat_file

    supplementarysecondary_reads_input=`bc <<< $(grep "secondary" $input_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $input_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    mapped_reads_input=`bc <<< $(grep "mapped (" $input_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$supplementarysecondary_reads_input`
    dup_reads_input=`grep "duplicates" $input_flagstat_file | sed -e 's/ + [[:digit:]]* duplicates$//'`
    dup_rate_input=$(echo "100*${dup_reads_input}/${mapped_reads_input}" | bc -l)

    if [[ -s $trimmomatic_table ]]
    then
      raw_reads_input=$(awk -v INPUT_NAME=$INPUT_NAME '{if ($1 == INPUT_NAME) print $0}' $trimmomatic_table | cut -f 2)
      trimmed_reads_input=`bc <<< $(grep "in total" $input_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$supplementarysecondary_reads_input` `grep "in total" $input_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//'`
      mapped_rate_input=$(echo "100*${mapped_reads_input}/${trimmed_reads_input}" | bc -l)
      trimmed_rate_input=$(echo "100*${trimmed_reads_input}/${raw_reads_input}" | bc -l)
    else
      raw_reads_input=`bc <<< $(grep "in total" $input_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$supplementarysecondary_reads_input`
      trimmed_reads_input="NULL"
      mapped_rate_input=$(echo "100*${mapped_reads_input}/${raw_reads_input}" | bc -l)
      trimmed_rate_input="NULL"
    fi
    # raw_reads_input=$(awk -v INPUT_NAME=$INPUT_NAME '{if ($1 == INPUT_NAME) print $0}' $trimmomatic_table | cut -f 2)
    # trimmed_reads_input=`grep "in total" $input_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//'`
    # mapped_reads_input=`grep "mapped (" $input_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
    # dup_reads_input=`grep "duplicates" $input_flagstat_file | sed -e 's/ + [[:digit:]]* duplicates$//'`
    # dup_rate_input=$(echo "100*${dup_reads_input}/${mapped_reads_input}" | bc -l)
    # mapped_rate_input=$(echo "100*${mapped_reads_input}/${trimmed_reads_input}" | bc -l)
    # trimmed_rate_input=$(echo "100*${trimmed_reads_input}/${raw_reads_input}" | bc -l)

    ## Finally, the number of singletons for paired-end data sets can be calculated using:
    singletons_input=`grep "singletons" $input_flagstat_file | sed -e 's/ + [[:digit:]]* singletons .*//'`

    # samtools view -b -F 3844 -q 5  ${INPUT_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}_INPUT.dedup.bam
    sambamba view -t $n -f bam -F "not unmapped and not secondary_alignment and not failed_quality_control and not duplicate and not supplementary and mapping_quality >= 5"  ${INPUT_BAM} > ${OUTPUT_DIR}/${SAMPLE_NAME}.Input.dedup.bam
    # samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}_INPUT.dedup.bam
    sambamba index -t $n ${OUTPUT_DIR}/${SAMPLE_NAME}.Input.dedup.bam
    # filtered_reads_input=`samtools flagstat  ${OUTPUT_DIR}/${SAMPLE_NAME}_INPUT.dedup.bam | grep "mapped (" | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
    filtered_reads_input=`sambamba flagstat -t $n ${OUTPUT_DIR}/${SAMPLE_NAME}.Input.dedup.bam | grep "mapped (" | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
    filtered_rate_input=$(echo "100*${filtered_reads_input}/${trimmed_reads_input}" | bc -l)
    # MT_reads_input=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE_NAME}_INPUT.dedup.bam MT)
    MT_reads_input=$(sambamba view -t $n -c ${OUTPUT_DIR}/${SAMPLE_NAME}.Input.dedup.bam MT)

    if [ -z $MT_reads_chip ] || [ $MT_reads_input -eq 0 ]
      then
        # MT_reads_input=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE_NAME}_INPUT.dedup.bam chrM)
        MT_reads_input=$(sambamba view -t $n -c ${OUTPUT_DIR}/${SAMPLE_NAME}.Input.dedup.bam chrM)
    fi
    MT_rate_input=$(echo "100*${MT_reads_input}/${filtered_reads_input}" | bc -l)
fi


# filtered_reads_chip=`samtools flagstat  ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam | grep "mapped (" | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
filtered_reads_chip=`sambamba flagstat -t $n $dedup_bam | grep "mapped (" | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
filtered_rate_chip=$(echo "100*${filtered_reads_chip}/${trimmed_reads_chip}" | bc -l)
# MT_reads_chip=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam MT)
MT_reads_chip=$(sambamba view -t $n -c $dedup_bam MT)

if [ -z $MT_reads_chip ] || [ $MT_reads_chip -eq 0 ]
  then
    # MT_reads_chip=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE_NAME}.dedup.bam chrM)
    MT_reads_chip=$(sambamba view -t $n -c $dedup_bam chrM)
fi
MT_rate_chip=$(echo "100*${MT_reads_chip}/${filtered_reads_chip}" | bc -l)

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
    plotFingerprint -b $dedup_bam $INPUT_BAM -bs ${bin_size} -l ${SAMPLE_NAME}.${CHIP_NAME} ${SAMPLE_NAME}.Input --JSDsample $INPUT_BAM --outQualityMetrics ${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.fingerprint.txt -plot ${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.fingerprint.png -p $n
    js_dist=`grep ${SAMPLE_NAME} ${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.fingerprint.txt | cut -f 8`
    chance_div=`grep ${SAMPLE_NAME} ${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.fingerprint.txt | cut -f 12`
fi

#4.     Calculating FRiP scores
nmb_peaks=$(wc -l ${CHIP_BED_FILE} | cut -f 1 -d " ")
reads_under_peaks=`samtools view -@ $n -c -L ${CHIP_BED_FILE} $dedup_bam`
frip=$(echo "${reads_under_peaks}/${filtered_reads_chip}" | bc -l)

#5. extract NSC and RSC from run_spp

nsc_chip=$(cut -f 9 ${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.crosscor | head -n 1)
rsc_chip=$(cut -f 10 ${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.crosscor | head -n 1)
quality_chip_num=$(cut -f 11 ${OUTPUT_DIR}/${SAMPLE_NAME}.${CHIP_NAME}.crosscor | head -n 1)

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

printf "genome_assembly\tChIP_type\ttreat_name\tctl_name\ttreat_raw_reads\ttreat_trimmed_reads\ttreat_trimmed_frac\ttreat_mapped_reads\ttreat_mapped_frac\ttreat_dup_reads\ttreat_dup_frac\ttreat_filtered_reads\ttreat_filtered_frac\ttreat_MT_reads\ttreat_MT_frac\tctl_raw_reads\tctl_trimmed_reads\tctl_trimmed_frac\tclt_mapped_reads\tctl_mapped_frac\tctl_dup_reads\tctl_dup_frac\tctl_filtered_reads\tctl_filtered_frac\tctl_MT_reads\tctl_Mt_frac\tnmb_peaks\treads_in_peaks\tfrip\ttreat_nsc\tctl_nsc\ttreat_rsc\tctl_rsc\ttreat_Quality\tctl_Quality\tsingletons\tjs_dist\tchance_div\n" > ${OUTPUT_DIR}/IHEC_metrics_chipseq_${SAMPLE_NAME}.${CHIP_NAME}.txt
LANG=C printf "%s\t%s\t%s\t%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\t%d\t%.4f\t%.4f\n" "${assembly}" "${CHIP_TYPE}" "${SAMPLE_NAME}.${CHIP_NAME}" "${INPUT_NAME}.Input" "$raw_reads_chip" "$trimmed_reads_chip" "$trimmed_rate_chip" "$mapped_reads_chip" "$mapped_rate_chip" "$dup_reads_chip" "$dup_rate_chip" "$filtered_reads_chip" "$filtered_rate_chip" "$MT_reads_chip" "$MT_rate_chip" "$raw_reads_input" "$trimmed_reads_input" "$trimmed_rate_input" "$mapped_reads_input" "$mapped_rate_input" "$dup_reads_input" "$dup_rate_input" "$filtered_reads_input" "$filtered_rate_input" "$MT_reads_input" "$MT_rate_input" "$nmb_peaks" "$reads_under_peaks" "$frip" "$nsc_chip" "$nsc_input"  "$rsc_chip" "$rsc_input" "$quality_chip" "$quality_input" "$singletons_chip" "$js_dist" "$chance_div">> ${OUTPUT_DIR}/IHEC_metrics_chipseq_${SAMPLE_NAME}.${CHIP_NAME}.txt
