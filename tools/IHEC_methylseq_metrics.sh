#!/usr/bin/env sh
set -e
set -o pipefail

# Get args
SAMPLE_NAME=$1
TABLE_OUTFILE=$2
TABLE_OUTFILE_ALL=$3
TARGET_FLAG=$4
COUNT=$5


if [ $TARGET_FLAG == 1 ]; then
  echo -e "sample\traw_reads\ttrimmed_reads\t%_survivalrate\taln_reads\tMappingEfficiency\tDuplicatedReads\t%aligned_duplicate\tDeduplicatedAlignRreads\t%_UsefulAlignRate\t%Proportion_Unique_filteredReads\tOntargetReads\t%Ontarget\t%onTargetvsRawRead\tGC_bias\tpUC19_meth\tlambdaConversion\thumanConversion\tmean_genomecoverage\tmedianCpGcoverage\t#_CG_1X\t#_CG_10X\t#_CG_30X" > $TABLE_OUTFILE
  if [ $COUNT -eq 0 ]; then 
    echo -e "sample\traw_reads\ttrimmed_reads\t%_survivalrate\taln_reads\tMappingEfficiency\tDuplicatedReads\t%aligned_duplicate\tDeduplicatedAlignRreads\t%_UsefulAlignRate\t%Proportion_Unique_filteredReads\tOntargetReads\t%Ontarget\t%onTargetvsRawRead\tGC_bias\tpUC19_meth\tlambdaConversion\thumanConversion\tmean_genomecoverage\tmedianCpGcoverage\t#_CG_1X\t#_CG_10X\t#_CG_30X" > $TABLE_OUTFILE_ALL
  fi
else
  echo -e "sample\traw_reads\ttrimmed_reads\t%_survivalrate\taln_reads\tMappingEfficiency\tDuplicatedReads\t%aligned_duplicate\tDeduplicatedAlignRreads\t%_UsefulAlignRate\t%Proportion_Unique_filteredReads\tGC_bias\tpUC19_meth\tlambdaConversion\thumanConversion\tmean_genomecoverage\tmedianCpGcoverage\t#_CG_1X\t#_CG_10X\t#_CG_30X" > $TABLE_OUTFILE
  if [ $COUNT -eq 0 ]; then 
    echo -e "sample\traw_reads\ttrimmed_reads\t%_survivalrate\taln_reads\tMappingEfficiency\tDuplicatedReads\t%aligned_duplicate\tDeduplicatedAlignRreads\t%_UsefulAlignRate\t%Proportion_Unique_filteredReads\tGC_bias\tpUC19_meth\tlambdaConversion\thumanConversion\tmean_genomecoverage\tmedianCpGcoverage\t#_CG_1X\t#_CG_10X\t#_CG_30X" > $TABLE_OUTFILE_ALL
  fi
fi




# Calculate the survival rate : #reads after trimming
rawReads=`cat trim/${SAMPLE_NAME}/*.trim.log |grep "Input"| sed s/"Input Read Pairs: "//g|sed s/"Both Surviving:"//g|awk '{sum+=$1;} END {printf "%d\n", sum}'`
trimmedReads=`cat trim/${SAMPLE_NAME}/*.trim.log |grep "Input"| sed s/"Input Read Pairs: "//g|sed s/"Both Surviving:"//g|awk '{sum+=$2;} END {printf "%d\n", sum}'`
a=`echo $trimmedReads` && b=`echo $rawReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && SurvivalRate=`echo $nr`;

# Mapping efficiency after trimming, IHEC required.
AlignedReads=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.readset_sorted.deduplication_report.txt |grep -e "Total number of alignments analysed"|awk '{print $NF}'`
a=`echo $AlignedReads` && b=`echo $trimmedReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && MappingEfficiency=`echo $nr`;

# De-duplicate rate after alignment
DuplicateReads=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.readset_sorted.deduplication_report.txt |grep -e "Total number duplicated alignments removed:"|awk -F "\t" '{print $2}' |awk '{print $1}'`
a=`echo $DuplicateReads` && b=`echo $AlignedReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && DuplicationRate=`echo $nr`;
a=`echo $DuplicateReads` && b=`echo $AlignedReads` && DeduplicatedAlignRreads=$(echo " ($b-$a)" |bc)  
a=`echo $DuplicateReads` && b=`echo $AlignedReads` && c=`echo $rawReads` && nr=$(echo "scale=4;( ($b-$a) / $c) * 100;" | bc) && UsefulAlignRate=`echo $nr`;

# Estimated average coverage, values for this metric should be above 12 for a single lane. IHEC required.
genomecoverage=`sed 1d alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.sorted.dedup.all.coverage.sample_summary | awk '{print $3}' | head -1`

# Specific for targeted capture methylome (MCC-Seq) data, obtain the on-target rate, the final useful proportion of reads over the raw reads.
if [ $TARGET_FLAG == 1 ]; then
  OntargetReads=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.sorted.dedup.ontarget.bam.flagstat|grep mapped|head -1|awk '{printf "%d",$1/2}'`
  a=`echo $OntargetReads` && b=`echo $DeduplicatedAlignRreads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && OntargetRate=`echo $nr`;
  a=`echo $OntargetReads` && b=`echo $rawReads` && nr=$(echo "scale=4;($a / $b) * 100;" | bc) && onTargetvsRawRead=`echo $nr`;
fi

# Proportion of uniquely aligned reads without duplicates and a quality score > 10; for reference epigenome, values for this metric should exceed 85%. IHEC required.
DedupMQ10filteredreads=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}.sorted.dedup.filtered_reads.counts.txt`
a=`echo $DedupMQ10filteredreads` && b=`echo $trimmedReads` && nr=$(echo "scale=4; (($a / $b) /2)* 100;" | bc) && DedupMQ10Proportion=`echo $nr`;

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
TotalmCHG=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}*/${SAMPLE_NAME}.*sorted_noRG_bismark_bt2_PE_report.txt|grep "Total methylated C's in CHG context:"|awk '{print $NF}'|awk '{sum+=$1;} END {printf "%d\n", sum}'`
TotalmCHH=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}*/${SAMPLE_NAME}.*sorted_noRG_bismark_bt2_PE_report.txt|grep "Total methylated C's in CHH context:"|awk '{print $NF}'|awk '{sum+=$1;} END {printf "%d\n", sum}'`
TotalumCHG=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}*/${SAMPLE_NAME}.*sorted_noRG_bismark_bt2_PE_report.txt|grep "Total unmethylated C's in CHG context:"|awk '{print $NF}'|awk '{sum+=$1;} END {printf "%d\n", sum}'`
TotalumCHH=`cat alignment/${SAMPLE_NAME}/${SAMPLE_NAME}*/${SAMPLE_NAME}.*sorted_noRG_bismark_bt2_PE_report.txt|grep "Total unmethylated C's in CHH context:"|awk '{print $NF}'|awk '{sum+=$1;} END {printf "%d\n", sum}'`
a=`echo $TotalmCHG` && b=`echo $TotalmCHH` && c=`echo $TotalumCHG` && d=`echo $TotalumCHH` &&  nr=$(echo "scale=4;( ($c+$d) / ($a+$b+$c+$d)) * 100;" | bc) && humanConversion=`echo $nr`;

if [ $TARGET_FLAG == 1 ]; then
  echo -e "${SAMPLE_NAME}\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$MappingEfficiency\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$DedupMQ10Proportion\t$OntargetReads\t$OntargetRate\t$onTargetvsRawRead\t$GCbias\t$puc19Meth\t$lambdaConversion\t$humanConversion\t$genomecoverage\t$MedianCpGcov\t$cg1x\t$cg10x\t$cg30x" >> $TABLE_OUTFILE
  echo -e "${SAMPLE_NAME}\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$MappingEfficiency\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$DedupMQ10Proportion\t$OntargetReads\t$OntargetRate\t$onTargetvsRawRead\t$GCbias\t$puc19Meth\t$lambdaConversion\t$humanConversion\t$genomecoverage\t$MedianCpGcov\t$cg1x\t$cg10x\t$cg30x" >> $TABLE_OUTFILE_ALL
else
  echo -e "${SAMPLE_NAME}\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$MappingEfficiency\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$DedupMQ10Proportion\t$GCbias\t$puc19Meth\t$lambdaConversion\t$humanConversion\t$genomecoverage\t$MedianCpGcov\t$cg1x\t$cg10x\t$cg30x" >> $TABLE_OUTFILE
  echo -e "${SAMPLE_NAME}\t$rawReads\t$trimmedReads\t$SurvivalRate\t$AlignedReads\t$MappingEfficiency\t$DuplicateReads\t$DuplicationRate\t$DeduplicatedAlignRreads\t$UsefulAlignRate\t$DedupMQ10Proportion\t$GCbias\t$puc19Meth\t$lambdaConversion\t$humanConversion\t$genomecoverage\t$MedianCpGcov\t$cg1x\t$cg10x\t$cg30x" >> $TABLE_OUTFILE_ALL
fi

