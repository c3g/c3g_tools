#!/usr/bin/env sh
set -euxo pipefail

# Get args
SAMPLE_LIST=$1
TABLE_OUTFILE=$2
TARGET_FLAG=$3

if [ $TARGET_FLAG == 1 ]; then
  echo -e "sample\traw_reads\ttrimmed_reads\t%_survivalrate\taligned_reads\t%_mapping_efficiency\tduplicated_reads\t%_duplication_rate\tdeduplicated_aligned_reads\t%_useful_aligned_rate\ton_target_reads\t%_on_target_rate\t%_on_target_vs_raw_reads\t%_lambda_conversion_rate\testimated_average_genome_coverage\t#_CG_1X\t#_CG_10X\t#_CG_30X" > $TABLE_OUTFILE
else
  echo -e "sample\traw_reads\ttrimmed_reads\t%_survivalrate\taligned_reads\t%_mapping_efficiency\tduplicated_reads\t%_duplication_rate\tdeduplicated_aligned_reads\t%_useful_aligned_rate\t%_lambda_conversion_rate\testimated_average_genome_coverage\t#_CG_1X\t#_CG_10X\t#_CG_30X" > $TABLE_OUTFILE
fi

for sample in `echo $SAMPLE_LIST | sed 's/,/ /g'`
do
  rawReads=`cat trim/$sample/*.trim.log |grep "Input"| sed s/"Input Read Pairs: "//g|sed s/"Both Surviving:"//g|awk '{sum+=$1;} END {printf "%d\n", sum*2}'`
  trimmedReads=`cat trim/$sample/*.trim.log |grep "Input"| sed s/"Input Read Pairs: "//g|sed s/"Both Surviving:"//g|awk '{sum+=$2;} END {printf "%d\n", sum*2}'`
  a=`echo $trimmedReads` && b=`echo $rawReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && SurvivalRate=`echo $nr`;

  # The number of aligned reads :
  samtools flagstat alignment/$sample/$sample.sorted.bam > alignment/$sample/$sample.sorted_flagstat.txt
  AlignedReads=`grep "mapped (" alignment/$sample/$sample.sorted_flagstat.txt | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
  a=`echo $AlignedReads` && b=`echo $trimmedReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && AlignedRate2=`echo $nr`;

  # The number of deduplicated aligned reads:
  samtools flagstat alignment/$sample/$sample.sorted.dedup.bam > alignment/$sample/$sample.sorted.dedup_flagstat.txt
  DeduplicatedAlignRreads=`grep "mapped (" alignment/$sample/$sample.sorted.dedup_flagstat.txt | sed -e 's/ + [[:digit:]]* mapped (.*)//'`

  DuplicateReads=$(echo " ($AlignedReads-$DeduplicatedAlignRreads)" | bc)
  DuplicationRate=$(echo "scale=4;(${DuplicateReads} / ${AlignedReads}) * 100" | bc -l)
  a=`echo $DuplicateReads` && b=`echo $AlignedReads` && c=`echo $rawReads` && nr=$(echo "scale=4;( ($b-$a) / $c) * 100;" | bc) && UsefulAlignRate=`echo $nr`;

  coverage=`sed 1d alignment/$sample/$sample.sorted.dedup.all.coverage.sample_summary | awk '{print $3}' | tail -1`

  # Check if the ontarget flagstat file exists (i.e. captured analysis or not)
  if [ $TARGET_FLAG == 1 ]; then
    OntargetReads=`cat alignment/$sample/$sample.sorted.dedup.ontarget.bam.flagstat|grep mapped|head -1|awk '{printf "%d",$1}'`
    a=`echo $OntargetReads` && b=`echo $DeduplicatedAlignRreads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && OntargetRate=`echo $nr`;
    a=`echo $OntargetReads` && b=`echo $rawReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && onTargetvsRawRead=`echo $nr`;
  fi

  lambdaConversion=`cat methylation_call/$sample/$sample*.profile.lambda.conversion.rate.tsv | awk '{print $1}'`

  #_CG_1X,#_CG_10X,#_CG_30X
  file=methylation_call/$sample/$sample*.profile.cgstats.txt
  cg1x=`cat  $file | awk -F"," '{print $1}'`
  cg10x=`cat $file | awk -F"," '{print $2}'`
  cg30x=`cat $file | awk -F"," '{print $3}'`

  if [ $TARGET_FLAG == 1 ]; then
    echo -e "$sample\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$AlignedRate2\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$OntargetReads\t$OntargetRate\t$onTargetvsRawRead\t$lambdaConversion\t$coverage\t$cg1x\t$cg10x\t$cg30x"
  else
    echo -e "$sample\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$AlignedRate2\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$lambdaConversion\t$coverage\t$cg1x\t$cg10x\t$cg30x"
  fi

done >> $TABLE_OUTFILE
