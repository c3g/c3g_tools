#!/usr/bin/env sh
set -euxo pipefail

# Get args
SAMPLE_NAME=$1
TABLE_OUTFILE=$2
TABLE_OUTFILE_ALL=$3
TARGET_FLAG=$4
COUNT=$5


if [ $TARGET_FLAG == 1 ]; then
  echo -e "sample\traw_reads\ttrimmed_reads\t%_survival_rate\taligned_reads\t%_mapping_efficiency\tduplicated_reads\t%_duplication_rate\tdeduplicated_aligned_reads\t%_useful_aligned_rate\t%_proportion_unique_filtered_reads_MAPQ>10\ton_target_dedup_reads\t%_on_target_rate\t%on_target_vs_raw_reads\tGC_bias\t%_pUC19_methylation_rate\t%_lambda_conversion_rate\t%_human_conversion\testimated_average_genome_coverage\tmedian_CpG_coverage\t#_CpG_1X\t#_CpG_10X\t#_CpG_30X\testimate_library_size\t#_onTrarget_CpG_1X\t#_onTrarget_CpG_10X\t#_onTrarget_CpG_15X\t#_onTrarget_CpG_20X\t#_onTrarget_CpG_25X\t#_onTrarget_CpG_30X\ton_target_raw_reads\t%_onTrarget_duplication_rate" > $TABLE_OUTFILE
  if [ $COUNT -eq 0 ]; then 
    echo -e "sample\traw_reads\ttrimmed_reads\t%_survival_rate\taligned_reads\t%_mapping_efficiency\tduplicated_reads\t%_duplication_rate\tdeduplicated_aligned_reads\t%_useful_aligned_rate\t%_proportion_unique_filtered_reads_MAPQ>10\ton_target_dedup_reads\t%_on_target_rate\t%on_target_vs_raw_reads\tGC_bias\t%_pUC19_methylation_rate\t%_lambda_conversion_rate\t%_human_conversion\testimated_average_genome_coverage\tmedian_CpG_coverage\t#_CpG_1X\t#_CpG_10X\t#_CpG_30X\testimate_library_size\t#_onTrarget_CpG_1X\t#_onTrarget_CpG_10X\t#_onTrarget_CpG_15X\t#_onTrarget_CpG_20X\t#_onTrarget_CpG_25X\t#_onTrarget_CpG_30X\ton_target_raw_reads\t%_onTrarget_duplication_rate" > $TABLE_OUTFILE_ALL
  fi
else
  echo -e "sample\traw_reads\ttrimmed_reads\t%_survival_rate\taligned_reads\t%_mapping_efficiency\tduplicated_reads\t%_duplication_rate\tdeduplicated_aligned_reads\t%_useful_aligned_rate\t%_proportion_unique_filtered_reads_MAPQ>10\tGC_bias\t%_pUC19_methylation_rate\t%_lambda_conversion_rate\t%_human_conversion\testimated_average_genome_coverage\tmedian_CpG_coverage\t#_CpG_1X\t#_CpG_10X\t#_CpG_30X\testimate_library_size" > $TABLE_OUTFILE
  if [ $COUNT -eq 0 ]; then 
    echo -e "sample\traw_reads\ttrimmed_reads\t%_survival_rate\taligned_reads\t%_mapping_efficiency\tduplicated_reads\t%_duplication_rate\tdeduplicated_aligned_reads\t%_useful_aligned_rate\t%_proportion_unique_filtered_reads_MAPQ>10\tGC_bias\t%_pUC19_methylation_rate\t%_lambda_conversion_rate\t%_human_conversion\testimated_average_genome_coverage\tmedian_CpG_coverage\t#_CpG_1X\t#_CpG_10X\t#_CpG_30X\testimate_library_size" > $TABLE_OUTFILE_ALL
  fi
fi

# Calculate the survival rate : #reads after trimming
rawReads=`cat trim/${SAMPLE_NAME}/*.trim.log |grep "Input"| sed s/"Input Read Pairs: "//g|sed s/"Both Surviving:"//g|awk '{sum+=$1;} END {printf "%d\n", sum*2}'`
trimmedReads=`cat trim/${SAMPLE_NAME}/*.trim.log |grep "Input"| sed s/"Input Read Pairs: "//g|sed s/"Both Surviving:"//g|awk '{sum+=$2;} END {printf "%d\n", sum*2}'`
a=`echo $trimmedReads` && b=`echo $rawReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && SurvivalRate=`echo $nr`;

# The number of aligned reads :
AlignedReads=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.raw.count | head -1 |awk '{printf "%d",$1}'`
a=`echo $AlignedReads` && b=`echo $trimmedReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && MappingEfficiency=`echo $nr`;

# The number of deduplicated aligned reads:
DeduplicatedAlignRreads=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.dedup.count | head -1 |awk '{printf "%d",$1}'`
DuplicateReads=$(echo " ($AlignedReads-$DeduplicatedAlignRreads)" | bc)
DuplicationRate=$(echo "scale=4;(${DuplicateReads} / ${AlignedReads}) * 100;" | bc -l)

a=`echo $DuplicateReads` && b=`echo $AlignedReads` && c=`echo $rawReads` && nr=$(echo "scale=4;( ($b-$a) / $c) * 100;" | bc) && UsefulAlignRate=`echo $nr`;

# Estimated average coverage, values for this metric should be above 12 for a single lane. IHEC required.
genomecoverage=`sed 1d alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.sorted.dedup.all.coverage.sample_summary | awk '{print $3}' | tail -1`

# Estimated library estimate library size 
est_lib_size=`grep -A1 ESTIMATED_LIBRARY_SIZE alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.sorted.dedup.metrics | tail -1 | cut -f10`

# Specific for targeted capture methylome (MCC-Seq) data, obtain the on-target rate (raw and dedup), the final useful proportion of reads over the raw reads.
if [ $TARGET_FLAG == 1 ]; then
  OntargetReadsDedup=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.onTarget.dedup.count | head -1 |awk '{printf "%d",$1}'`
  a=`echo $OntargetReadsDedup` && b=`echo $DeduplicatedAlignRreads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && OntargetRate=`echo $nr`;
  a=`echo $OntargetReadsDedup` && b=`echo $rawReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && onTargetvsRawRead=`echo $nr`;
  OntargetReadsRaw=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.onTarget.raw.count | head -1 |awk '{printf "%d",$1}'`
  a=`echo $OntargetReadsDedup` && b=`echo $OntargetReadsRaw` && nr=$(echo "scale=4;1-($a / $b) * 100;" | bc) && onTargetDedupRate=`echo $nr`;
  file=methylation_call/${SAMPLE_NAME}/${SAMPLE_NAME}.readset_sorted.dedup.CpG_profile.strand.combined.on_target.count
  ot_cg1x=`tail -1  $file | awk '{print $2}'`
  ot_cg10x=`tail -1 $file | awk '{print $3}'`
  ot_cg15x=`tail -1 $file | awk '{print $4}'`
  ot_cg20x=`tail -1 $file | awk '{print $5}'`
  ot_cg25x=`tail -1 $file | awk '{print $6}'`
  ot_cg30x=`tail -1 $file | awk '{print $7}'`
  
fi

# Proportion of uniquely aligned reads without duplicates and a quality score > 10; for reference epigenome, values for this metric should exceed 85%. IHEC required.
DedupMQ10filteredreads=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.sorted.dedup.filtered_reads.counts.txt`
a=`echo $DedupMQ10filteredreads` && b=`echo $trimmedReads` && nr=$(echo "scale=4; (($a / $b) )* 100;" | bc) && DedupMQ10Proportion=`echo $nr`;

# Median CpG coverage, values for this metrics should exceed 2. IHEC required.
MedianCpGcov=`cat methylation_call/${SAMPLE_NAME}/${SAMPLE_NAME}.readset_sorted.dedup.median_CpG_coverage.txt|head -1 |awk '{print $NF}'`

# GC bias, the absolute value of the correlation score should not exceed 0.4. IHEC required.
GCbias=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.sorted.dedup.GCBias_all.txt|tail -1`

# Report for the lambda conversion rate, estimate the bisulfite conversion rate, estimated from lambda genome reads; Values for this metrics should exceed 97% in both human and lambda. IHEC required.
lambdaConversion=`cat methylation_call/${SAMPLE_NAME}/${SAMPLE_NAME}.*profile.lambda.conversion.rate.tsv | awk '{print $1}'`

# Calculate for the puc19 methylation level, estimate of the over-bisulfite conversion rate.
puc19Meth=`cat methylation_call/${SAMPLE_NAME}/${SAMPLE_NAME}.*pUC19*.txt`

#_CG_1X,#_CG_10X,#_CG_30X
file=methylation_call/${SAMPLE_NAME}/${SAMPLE_NAME}.*profile.cgstats.txt
cg1x=`cat  $file | awk -F"," '{print $1}'`
cg10x=`cat $file | awk -F"," '{print $2}'`
cg30x=`cat $file | awk -F"," '{print $3}'`

# Calculate for the human conversion rate, estimated from bismark alignment reports, IHEC required.
TotalmCHG=$(for report in `find alignment/${SAMPLE_NAME} -name "*sorted_noRG_bismark_bt2_PE_report.txt"`; do cat $report | grep "Total methylated C's in CHG context:" | awk '{print $NF}'; done | awk '{sum+=$1;} END {printf "%d\n", sum}')
TotalmCHH=$(for report in `find alignment/${SAMPLE_NAME} -name "*sorted_noRG_bismark_bt2_PE_report.txt"`; do cat $report | grep "Total methylated C's in CHH context:" | awk '{print $NF}'; done | awk '{sum+=$1;} END {printf "%d\n", sum}')
TotalumCHG=$(for report in `find alignment/${SAMPLE_NAME} -name "*sorted_noRG_bismark_bt2_PE_report.txt"`; do cat $report | grep "Total unmethylated C's in CHG context:" | awk '{print $NF}'; done | awk '{sum+=$1;} END {printf "%d\n", sum}')
TotalumCHH=$(for report in `find alignment/${SAMPLE_NAME} -name "*sorted_noRG_bismark_bt2_PE_report.txt"`; do cat $report | grep "Total unmethylated C's in CHH context:" | awk '{print $NF}'; done | awk '{sum+=$1;} END {printf "%d\n", sum}')
a=`echo $TotalmCHG` && b=`echo $TotalmCHH` && c=`echo $TotalumCHG` && d=`echo $TotalumCHH` &&  nr=$(echo "scale=4;( ($c+$d) / ($a+$b+$c+$d)) * 100;" | bc) && humanConversion=`echo $nr`;

if [ $TARGET_FLAG == 1 ]; then
  echo -e "${SAMPLE_NAME}\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$MappingEfficiency\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$DedupMQ10Proportion\t$OntargetReadsDedup\t$OntargetRate\t$onTargetvsRawRead\t$GCbias\t$puc19Meth\t$lambdaConversion\t$humanConversion\t$genomecoverage\t$MedianCpGcov\t$cg1x\t$cg10x\t$cg30x\t$est_lib_size\t$ot_cg1x\t$ot_cg10x\t$ot_cg15x\t$ot_cg20x\t$ot_cg25x\t$ot_cg30x\t$OntargetReadsRaw\t$onTargetDedupRate" >> $TABLE_OUTFILE
  echo -e "${SAMPLE_NAME}\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$MappingEfficiency\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$DedupMQ10Proportion\t$OntargetReadsDedup\t$OntargetRate\t$onTargetvsRawRead\t$GCbias\t$puc19Meth\t$lambdaConversion\t$humanConversion\t$genomecoverage\t$MedianCpGcov\t$cg1x\t$cg10x\t$cg30x\t$est_lib_size\t$ot_cg1x\t$ot_cg10x\t$ot_cg15x\t$ot_cg20x\t$ot_cg25x\t$ot_cg30x\t$OntargetReadsRaw\t$onTargetDedupRate" >> $TABLE_OUTFILE_ALL
else
  echo -e "${SAMPLE_NAME}\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$MappingEfficiency\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$DedupMQ10Proportion\t$GCbias\t$puc19Meth\t$lambdaConversion\t$humanConversion\t$genomecoverage\t$MedianCpGcov\t$cg1x\t$cg10x\t$cg30x\t$est_lib_size" >> $TABLE_OUTFILE
  echo -e "${SAMPLE_NAME}\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$MappingEfficiency\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$DedupMQ10Proportion\t$GCbias\t$puc19Meth\t$lambdaConversion\t$humanConversion\t$genomecoverage\t$MedianCpGcov\t$cg1x\t$cg10x\t$cg30x\t$est_lib_size" >> $TABLE_OUTFILE_ALL
fi

