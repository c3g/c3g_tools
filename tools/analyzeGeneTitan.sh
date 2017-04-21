#!/usr/bin/env sh

# Exit immediately on error
set -eu -o pipefail
###bash script to test the SOP of the Axiom gene titan following the affy best practices
### mathieu bourgey - 2016/08/26 - mathieu.bourgey@mcgill.ca
VERSION=1.0

function usage {

            
            echo "Usage: `basename $0` [option] <value>"
            echo ""
            echo "analyzeGT.sh - launches the pipeline to process Gene titan data"
            echo ""
            echo "Options:"
            echo "####modules"
            echo "-A                    APT module (default mugqic_dev/AffymetricxApt/1.18.0)"
            echo "-M                    MUGQIC_TOOLS module (default mugqic_dev/mugqic_tools/2.1.6-beta)"
            echo "-R                    R module (default mugqic_dev/R_Bioconductor/3.2.3_3.2)"
            echo "-T                    thread the processing using 10 parallel thread for long steps (need multicpu available)"
            echo ""
            echo "####Path and files"
            echo "-c                    cell files folder (Mandatory)"
            echo "-l                    Chip library folder (containing xml files) (Mandatory)"
            echo "-o                    output folder (Mandatory)"
            echo "-m                    master list file containing the full list of cell files (Mandatory)"
            echo ""
            echo "####Generic parameters"
            echo "-a                    Axiom Array name (Mandatory)"
            echo "-r                    Axiom Array  revision name (Mandatory)"
            echo "-d                    DQC threshold for filtering [0-1] (default 0.82)"
            echo "-s                    Sample call rate threshold for filtering [0-100] (default 97)"
            echo "-p                    Plate pass rate threshold for filtering [0-100] (default 95)"
            echo "-C                    Average call rate for passing samples per batch [0-100] (default 98.5)" 
            echo "-S                    Specie of the sample (default human)"
            echo "-n                    number of SNP to plot for each category (default 6)"
            echo ""
            echo "#####Ps-classification  specific parameters (for experts only)"
            echo "-b                    Threshold for call rate. Default is 95.0"
            echo "-e                    Threshold for FLD.  Default value: 3.6"
            echo "-f                    Threshold for HetSO. Default value: -0.1"
            echo "-g                    Threshold for HetSO OTV. Default value: '-0.3"
            echo "-i                    Threshold for HomRO when SNP has 1 genotype. Default value: 0.6"
            echo "-j                    Threshold for HomRO when SNP has 2 genotype. Default value: 0.3"
            echo "-k                    Threshold for HomRO when SNP has 3 genotype. Default value: -0.9"
            echo "-l                    True if HomRO metric to be used in classification."
            echo "-n                    True if Hom Het metric to be used in classification."
            echo "-q                    Threshold for number of minor alleles. Default value: 2.0"
            echo "-u                    Priority order of probeset conversion types when performing probeset selection. Default value: PolyHighResolution,NoMinorHom,OTV,MonoHighResolution,CallRateBelowThreshold"
            echo "-v                    List of categories whose SNPs will be output as recommended. Default value is PolyHighResolution,MonoHighResolution,NoMinorHom,Hemizygous"
exit 0

}


##Default args
APT_MODULE=mugqic/AffymetrixApt/1.18.0
MUGIC_TOOLS_MODULE=mugqic/mugqic_tools/2.1.8
#r module should have the SNPolisher package installed
R_MODULE=mugqic_dev/R_Bioconductor/3.2.3_3.2
DQC_THRESHOLD=0.82
CALL_RATE_THRESHOLD=97
PLATE_PASS_RATE_THRESHOLD=95
AVERAGE_PLATE_CALL_RATE=98.5
SPECIE=human
OUTPUT_SNP_NUMBER=6
THREADED=0
CR_CUTOFF=95.00
FLD_CUTOFF=3.6
HET_SO_CUTOFF=-0.1
HET_SO_OTV_CUTOFF=-0.3
HOM_RO_1_CUTOFF=0.6
HOM_RO_2_CUTOFF=0.3
HOM_RO_3_CUTOFF=0.9
HOM_RO_BOOL=True
HOM_HET_BOOL=True
NMA_CUTOFF=2.0
PRIORITY=PolyHighResolution,NoMinorHom,OTV,MonoHighResolution,CallRateBelowThreshold
RECOMMENDED=PolyHighResolution,MonoHighResolution,NoMinorHom,Hemizygous


if [ $# -lt 12 ]
then
	usage
	exit 0
fi

##get arguments
while getopts "A:M:R:c:l:o:m:a:r:d:s:p:C:S:n:b:e:f:g:i:j:k:l:n:q:u:v:Th" OPT
do
    case "$OPT" in
        A) 
           APT_MODULE=$OPTARG
           echo "set-up APT module as: $APT_MODULE"
           ;;
	M) 
           MUGIC_TOOLS_MODULE=$OPTARG
           echo "set-up mugqic_tools module as: $MUGIC_TOOLS_MODULE"
           ;;
        R)
           R_MODULE=$OPTARG
           echo "set-up R module as: $R_MODULE"
           ;;
        c)
           CEL_PATH=$OPTARG
           echo "set-up cell path as: $CEL_PATH"
           ;;
        l)
           ANALYSIS_FILES_DIR=$OPTARG
           echo "set-up analysis file directory as: $ANALYSIS_FILES_DIR"
           ;;
        o)
           OUTDIR=$OPTARG
           echo "set-up output directory as: $OUTDIR"
           ;;
        m)
           MASTER_LIST=$OPTARG
           echo "set-up master list file as: $MASTER_LIST"
           ;;
        a)
           AXIOM_ARRAY_NAME=$OPTARG
           echo "set-up Axiom array name as: $AXIOM_ARRAY_NAME"
           ;;
        r)
           AXIOM_ARRAY_REV=$OPTARG
           echo "set-up Axiom array revision version as: $AXIOM_ARRAY_REV"
           ;;
        d)
           DQC_THRESHOLD=$OPTARG
           echo "set-up DishQC as: $DQC_THRESHOLD"
           ;;
        s)
          CALL_RATE_THRESHOLD=$OPTARG
          echo "set-up call rate threshold as: $CALL_RATE_THRESHOLD"
           ;;
        p)
           PLATE_PASS_RATE_THRESHOLD=$OPTARG
           echo "set-up Plate pass rate threshold as: $PLATE_PASS_RATE_THRESHOLD"
           ;;
        C)
           AVERAGE_PLATE_CALL_RATE=$OPTARG
           echo "set-up Average Plate Call Rate Threshold as: $AVERAGE_PLATE_CALL_RATE"
           ;;
        S) 
           SPECIE=$OPTARG
           echo "set-up specie as: $PECIE"
           ;;
        n)
           OUTPUT_SNP_NUMBER=$OPTARG
           echo "set-upNumber of SNP in output as: $OUTPUT_SNP_NUMBER"
           ;;
	T)
           THREADED=1
           echo "set-up multithread mode as: on (10 threads) "
           ;;
        b)
           CR_CUTOFF=$OPTARG
           ;;
        e)
           FLD_CUTOFF=$OPTARG
           ;;
	f)
           HET_SO_CUTOFF=$OPTARG
           ;;
	g)
           HET_SO_OTV_CUTOFF=$OPTARG
           ;;
        i)
           HOM_RO_1_CUTOFF=$OPTARG
           ;;
	j)
           HOM_RO_2_CUTOFF=$OPTARG
           ;;
        k)
           HOM_RO_3_CUTOFF=$OPTARG
           ;;
        l)
           HOM_RO_BOOL=$OPTARG
           ;;
	n)
           HOM_HET_BOOL=$OPTARG
           ;;
	q)
           NMA_CUTOFF=$OPTARG
           ;;
        u)
           PRIORITY=$OPTARG
           ;;
	v)
           RECOMMENDED=$OPTARG
           ;;
	h)
           usage
           exit
           ;;
        
    esac

done


REF_ID=$(date +"%F_%H-%M-%S")

#load APT 
module load ${APT_MODULE} ${MUGIC_TOOLS_MODULE} ${R_MODULE}

## dev_argument
R_TOOLS=~/work/repo/mugqic_tools/R-tools/

echo "$THREADED"
if [[ $THREADED == 1 ]]
then
echo "Running multithreaded"
fi

mkdir -p ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/  ${OUTDIR}/GENO_QC_${REF_ID}/

#######
#Step 0
#output standard information in sumary files
#######
BATCH_NAME=$(basename CEL_PATH)
echo "##########################" >  ${OUTDIR}/summary_${REF_ID}.txt
echo "Analysis summary" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Batch Name: ${BATCH_NAME}"  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Array Pacakage Name: ${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Array Type Name: ${AXIOM_ARRAY_NAME}" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Array Display Name: NA"  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Workflow Type: Best Practices Workflow"  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Date Created: ${REF_ID}"  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "" >>  ${OUTDIR}/summary_${REF_ID}.txt




########
## step1 Group Samples into Batches
## create file containing the list of cel files path
########
echo "-----------------------"
echo "Runing Step 1 - Sample Grouping"
echo "-----------------------"

echo cel_files > ${OUTDIR}/cel_list1_${REF_ID}.txt
ls ${CEL_PATH}/*.CEL >>  ${OUTDIR}/cel_list1_${REF_ID}.txt

continue1=$(grep -v "^cel_files$" ${OUTDIR}/cel_list1_${REF_ID}.txt | wc -l | cut -d\  -f1)

###sample summary
echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Sample summary" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Number of imput sample: ${continue1}"  >>  ${OUTDIR}/summary_${REF_ID}.txt

if [[ $continue1 == 0 ]]
then

echo "No sample founds"
exit 0

fi

########
## step2 Generate the Sample “DQC” Values Using APT
## DQC values are produced by the program apt-geno-qc
########
echo "-----------------------"
echo "Runing Step 2 - DQC estimation"
echo "-----------------------"
echo "CMD: apt-geno-qc --analysis-files-path ${ANALYSIS_FILES_DIR} \
--xml-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.apt-geno-qc.AxiomQC1.xml \
--cel-files ${OUTDIR}/cel_list1_${REF_ID}.txt \
--out-file ${OUTDIR}/GENO_QC_${REF_ID}/apt-geno-qc_${REF_ID}.txt \
--log-file ${OUTDIR}/GENO_QC_${REF_ID}/apt-geno-qc_${REF_ID}.log"

apt-geno-qc --analysis-files-path ${ANALYSIS_FILES_DIR} \
--xml-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.apt-geno-qc.AxiomQC1.xml \
--cel-files ${OUTDIR}/cel_list1_${REF_ID}.txt \
--out-file ${OUTDIR}/GENO_QC_${REF_ID}/apt-geno-qc_${REF_ID}.txt \
--log-file ${OUTDIR}/GENO_QC_${REF_ID}/apt-geno-qc_${REF_ID}.log

########
## Step 3: Conduct Sample QC on DQC
## Remove samples with a DQC value less than a given DQC threshold (default of 0.82)
####
echo "-----------------------"
echo "Runing Step 3 - DQC filtering"
echo "-----------------------"
echo "CMD: Rscript ${R_TOOLS}/filterAxiom.R -s genoqc \
-q ${OUTDIR}/GENO_QC_${REF_ID}/apt-geno-qc_${REF_ID}.txt \
-c ${OUTDIR}/cel_list1_${REF_ID}.txt \
-d ${DQC_THRESHOLD} \
-o ${OUTDIR}/cel_list2_${REF_ID}.txt"

Rscript ${R_TOOLS}/filterAxiom.R -s genoqc \
-q ${OUTDIR}/GENO_QC_${REF_ID}/apt-geno-qc_${REF_ID}.txt \
-c ${OUTDIR}/cel_list1_${REF_ID}.txt \
-d ${DQC_THRESHOLD} \
-o ${OUTDIR}/cel_list2_${REF_ID}.txt

continue2=$(grep -v "^cel_files$" ${OUTDIR}/cel_list2_${REF_ID}.txt | wc -l | cut -d\  -f1)

###sample summary
echo "Sample passing DQC: ${continue2} out of ${continue1}"  >>  ${OUTDIR}/summary_${REF_ID}.txt


if [[ $continue2 == 0 ]]
then

echo "No sample remains after DQC filtering"
exit 0

fi

########
## Step 4: Generate Sample QC Call Rates
########
echo "-----------------------"
echo "Runing Step 4 - Call Rates QC"
echo "-----------------------"

if [[ $THREADED == 1 ]]
then

echo "multithreading 10cpu"

for i in `seq 1 10` 
do 

echo "CMD: apt-genotype-axiom --log-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/apt-genotype-axiom.step4_${REF_ID}.log \
--arg-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}_96orMore_Step1.${AXIOM_ARRAY_REV}.apt-probeset-genotype.AxiomGT1.xml \
--analysis-files-path ${ANALYSIS_FILES_DIR} \
--out-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/tmp${i}/ \
--probeset-ids ${ANALYSIS_FILES_DIR}/chunked_10_ps/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.step1.part${i}.ps \
--cel-files ${OUTDIR}/cel_list2_${REF_ID}.txt"


mkdir -p ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/tmp${i}
apt-genotype-axiom --log-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/apt-genotype-axiom.step4_${REF_ID}.log \
--arg-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}_96orMore_Step1.${AXIOM_ARRAY_REV}.apt-probeset-genotype.AxiomGT1.xml \
--analysis-files-path ${ANALYSIS_FILES_DIR} \
--out-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/tmp${i}/ \
--probeset-ids ${ANALYSIS_FILES_DIR}/chunked_10_ps/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.step1.part${i}.ps \
--cel-files ${OUTDIR}/cel_list2_${REF_ID}.txt & pids+=($!)


done

wait "${pids[@]}" 


echo "CMD: Rscript ${R_TOOLS}/filterAxiom.R -s merge \
-n 10 \
-o ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}"


Rscript ${R_TOOLS}/filterAxiom.R -s merge \
-n 10 \
-o ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}



else

echo "CMD: apt-genotype-axiom --log-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/apt-genotype-axiom.step4_${REF_ID}.log \
--arg-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}_96orMore_Step1.${AXIOM_ARRAY_REV}.apt-probeset-genotype.AxiomGT1.xml \
--analysis-files-path ${ANALYSIS_FILES_DIR} \
--out-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/ \
--cel-files ${OUTDIR}/cel_list2_${REF_ID}.txt"

apt-genotype-axiom --log-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/apt-genotype-axiom.step4_${REF_ID}.log \
--arg-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}_96orMore_Step1.${AXIOM_ARRAY_REV}.apt-probeset-genotype.AxiomGT1.xml \
--analysis-files-path ${ANALYSIS_FILES_DIR} \
--out-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/ \
--cel-files ${OUTDIR}/cel_list2_${REF_ID}.txt

fi

########
## Step 5: QC the Samples Based on QC Call Rate
## Remove samples with a QC call rate value less than a given threshold (default of 97%)
#######
echo "-----------------------"
echo "Runing Step 5 - Sample QC"
echo "-----------------------"
echo "CMD: Rscript ${R_TOOLS}/filterAxiom.R -s sampleqc \
-q ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.report.txt \
-c ${OUTDIR}/cel_list2_${REF_ID}.txt \
-d ${CALL_RATE_THRESHOLD} \
-o ${OUTDIR}/cel_list3_${REF_ID}.txt"


Rscript ${R_TOOLS}/filterAxiom.R -s sampleqc \
-q ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.report.txt \
-c ${OUTDIR}/cel_list2_${REF_ID}.txt \
-d ${CALL_RATE_THRESHOLD} \
-o ${OUTDIR}/cel_list3_${REF_ID}.txt

continue5=$(grep -v "^cel_files$" ${OUTDIR}/cel_list3_${REF_ID}.txt | wc -l | cut -d\  -f1)

###sample summary
echo "Sample passing DQC and QC CR: ${continue5} out of ${continue1}"  >>  ${OUTDIR}/summary_${REF_ID}.txt


if [[ $continue5 == 0 ]]
then

echo "No sample remains after Sample QC"
exit 0

fi


########
## Step 6: QC the Plates
## remove a plate if average call rate of passing samples value less than a given threshold (default of 98.5)
#######
echo "-----------------------"
echo "Runing Step 6 - Plates QC"
echo "-----------------------"
echo "CMD: Rscript ${R_TOOLS}/filterAxiom.R -s plateqc \
-q ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.report.txt \
-t ${OUTDIR}/GENO_QC_${REF_ID}/apt-geno-qc_${REF_ID}.txt \
-c ${OUTDIR}/cel_list3_${REF_ID}.txt \
-i ${OUTDIR}/cel_list2_${REF_ID}.txt \
-d ${PLATE_PASS_RATE_THRESHOLD} \
-a ${AVERAGE_PLATE_CALL_RATE} \
-p ${MASTER_LIST} \
-o ${OUTDIR}/cel_list4_${REF_ID}.txt"

Rscript ${R_TOOLS}/filterAxiom.R -s plateqc \
-q ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.report.txt \
-t ${OUTDIR}/GENO_QC_${REF_ID}/apt-geno-qc_${REF_ID}.txt \
-c ${OUTDIR}/cel_list3_${REF_ID}.txt \
-i ${OUTDIR}/cel_list2_${REF_ID}.txt \
-d ${PLATE_PASS_RATE_THRESHOLD} \
-a ${AVERAGE_PLATE_CALL_RATE} \
-p ${MASTER_LIST} \
-o ${OUTDIR}/cel_list4_${REF_ID}.txt


continue6=$(grep -v "^cel_files$"  ${OUTDIR}/cel_list4_${REF_ID}.txt | wc -l | cut -d\  -f1)
perc6=$(echo "scale=5; (${continue6} / ${continue1})*100" | bc -l)
sampleNum=$(echo "${continue1} - ${continue6}" | bc -l)

###sample summary
echo "Sample passing DQC and QC CR and Plate QC: ${continue6} out of ${continue1} (${perc6}%)"  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Sample did not pass: ${sampleNum}"  >>  ${OUTDIR}/summary_${REF_ID}.txt


if [[ $continue6 == 0 ]]
then

echo "No sample remains after Plate QC"
exit 0

fi


########
## Step 7: Genotype Passing Sample
########
echo "-----------------------"
echo "Runing Step 7 - Genotyping"
echo "-----------------------"
if [[ $THREADED == 1 ]]
then

echo "multithreading 10cpu"

for i in `seq 1 10` 
do 

echo "CMD: apt-genotype-axiom --log-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/apt-genotype-axiom.step7_${REF_ID}.${i}.log \
--arg-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}_96orMore_Step2.${AXIOM_ARRAY_REV}.apt-axiom-genotype.Biallelic.AxiomGT1.apt2.xml \
--analysis-files-path ${ANALYSIS_FILES_DIR} \
--out-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/tmp${i}/ \
--summaries \
--write-models \
--cc-chp-output \
--dual-channel-normalization TRUE \
--probeset-ids ${ANALYSIS_FILES_DIR}/chunked_10_ps/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.Biallelic.part${i}.ps \
--cel-files ${OUTDIR}/cel_list4_${REF_ID}.txt"

mkdir -p ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/tmp${i}
apt-genotype-axiom --log-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/apt-genotype-axiom.step7_${REF_ID}.log \
--arg-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}_96orMore_Step2.${AXIOM_ARRAY_REV}.apt-axiom-genotype.Biallelic.AxiomGT1.apt2.xml \
--analysis-files-path ${ANALYSIS_FILES_DIR} \
--out-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/tmp${i}/ \
--summaries \
--write-models \
--cc-chp-output \
--dual-channel-normalization TRUE \
--probeset-ids ${ANALYSIS_FILES_DIR}/chunked_10_ps/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.Biallelic.part${i}.ps \
--cel-files ${OUTDIR}/cel_list4_${REF_ID}.txt & pids+=($!)


done

wait "${pids[@]}" 

echo "CMD: Rscript ${R_TOOLS}/filterAxiom.R -s merge \
-n 10 \
-o ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}"


Rscript ${R_TOOLS}/filterAxiom.R -s merge \
-n 10 \
-o ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}


else

echo "CMD: apt-genotype-axiom --log-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/apt-genotype-axiom.step7_${REF_ID}.log \
--arg-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}_96orMore_Step2.${AXIOM_ARRAY_REV}.apt-axiom-genotype.Biallelic.AxiomGT1.apt2.xml \
--analysis-files-path ${ANALYSIS_FILES_DIR} \
--out-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/ \
--summaries \
--write-models \
--cc-chp-output \
--dual-channel-normalization TRUE \
--cel-files ${OUTDIR}/cel_list4_${REF_ID}.txt"

apt-genotype-axiom --log-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/apt-genotype-axiom.step7_${REF_ID}.log \
--arg-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}_96orMore_Step2.${AXIOM_ARRAY_REV}.apt-axiom-genotype.Biallelic.AxiomGT1.apt2.xml \
--analysis-files-path ${ANALYSIS_FILES_DIR} \
--out-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/ \
--summaries \
--write-models \
--cc-chp-output \
--dual-channel-normalization TRUE \
--cel-files ${OUTDIR}/cel_list4_${REF_ID}.txt

fi

continue7=$(grep "^probeset_id" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.calls.txt | tr '\t' '\n' | grep -v  "probeset_id"  | wc -l)
echo "Number of Samples genotyped: ${continue7}"  >>  ${OUTDIR}/summary_${REF_ID}.txt
average_call_rate=$(grep -v "^#" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.report.txt | awk ' BEGIN {cr=0;ln=0} NR > 1 {cr+=$3;ln++} END { print cr/ln} ')
)
echo "Average QC CR for the passing samples: ${average_call_rate}"  >>  ${OUTDIR}/summary_${REF_ID}.txt
gender=$(grep -v "^#" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.report.txt | awk ' BEGIN {ma=0;fe=0;unk=0} NR > 1 {if ($2 == "male") {ma++} else if ($2 == "female"){fe++} else {unk++}} END { print "male=" ma " female=" fe " unknown=" unk} ')
echo "Gender Calls Counts: ${gender}"  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo ""  >>  ${OUTDIR}/summary_${REF_ID}.txt

########
## Step 8A: Run Ps-Metrics
########
echo "-----------------------"
echo "Runing Step 8A - SNP metrics"
echo "-----------------------"
echo "CMD: ps-metrics --posterior-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.snp-posteriors.txt \
--call-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.calls.txt \
--metrics-file  ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/PS_metrics.txt"


ps-metrics --posterior-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.snp-posteriors.txt \
--call-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.calls.txt \
--metrics-file  ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/PS_metrics.txt

########
## Step 8B: Run Ps_Classification
########
echo "-----------------------"
echo "Runing Step 8B - SNP Classification"
echo "-----------------------"
echo "CMD: ps-classification --species-type ${SPECIE} \
--metrics-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/PS_metrics.txt \
--output-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/ \
--ps2snp-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.ps2snp_map.ps \
--recommended $RECOMMENDED \
--cr-cutoff $CR_CUTOFF \
--fld-cutoff $FLD_CUTOFF \
--het-so-cutoff $HET_SO_CUTOFF \
--het-so-otv-cutoff $HET_SO_OTV_CUTOFF \
--hom-ro-1-cutoff $HOM_RO_1_CUTOFF \
--hom-ro-2-cutoff $HOM_RO_2_CUTOFF \
--hom-ro-3-cutoff $HOM_RO_3_CUTOFF \
--hom-ro $HOM_RO_BOOL \
--hom-het $HOM_HET_BOOL \
--num-minor-allele-cutoff $NMA_CUTOFF \
--priority-order $PRIORITY"


ps-classification --species-type ${SPECIE} \
--metrics-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/PS_metrics.txt \
--output-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/ \
--ps2snp-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.ps2snp_map.ps \
--recommended $RECOMMENDED \
--cr-cutoff $CR_CUTOFF \
--fld-cutoff $FLD_CUTOFF \
--het-so-cutoff $HET_SO_CUTOFF \
--het-so-otv-cutoff $HET_SO_OTV_CUTOFF \
--hom-ro-1-cutoff $HOM_RO_1_CUTOFF \
--hom-ro-2-cutoff $HOM_RO_2_CUTOFF \
--hom-ro-3-cutoff $HOM_RO_3_CUTOFF \
--hom-ro $HOM_RO_BOOL \
--hom-het $HOM_HET_BOOL \
--num-minor-allele-cutoff $NMA_CUTOFF \
--priority-order $PRIORITY


########
## Step 8C: Run Off Target Variant (OTV) caller
########
echo "-----------------------"
echo "Runing Step 8C - Off Target Variant (OTV) caller"
echo "-----------------------"
echo "CMD: otv-caller --summary-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.summary.txt \
--posterior-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.snp-posteriors.txt \
--call-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.calls.txt  \
--confidence-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.confidences.txt \
--pid-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/OffTargetVariant.ps \
--output-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/OTV/"



otv-caller --summary-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.summary.txt \
--posterior-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.snp-posteriors.txt \
--call-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.calls.txt  \
--confidence-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.confidences.txt \
--pid-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/OffTargetVariant.ps \
--output-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/OTV/



########
## Step 8D: RunSNP_polisher to generate plots
########
echo "-----------------------"
echo "Runing Step 8D - Generating plots with SNP_polisher "
echo "-----------------------"
echo "CMD: Rscript ${R_TOOLS}/filterAxiom.R -s plotSNP \
-d ${OUTPUT_SNP_NUMBER} \
-o ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}"


Rscript ${R_TOOLS}/filterAxiom.R -s plotSNP \
-d ${OUTPUT_SNP_NUMBER} \
-o ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}


echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "SNP Metrics Summary" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
phr=$(grep -c -v "^probeset_id" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/PolyHighResolution.ps))
nmh=$(grep -c -v "^probeset_id" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/NoMinorHom.ps)
mhr=$(grep -c -v "^probeset_id" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/MonoHighResolution.ps)
other=$(grep -c -v "^probeset_id" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/Other.ps)
crbt=$(grep -c -v "^probeset_id" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/CallRateBelowThreshold.ps)
hem=$(grep -c -v "^probeset_id" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/Hemizygous.ps)
otv=$(grep -c -v "^probeset_id" ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/OffTargetVariant.ps)
tot=$(echo "$phr + $nmh + $mhr + $other + $crbt + $hem + $otv" | bc -l)
echo "Number of SNPs: ${tot}" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo ""  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo -e "ConversionType\tCount\tPercentage"
perPHR=$(echo "scale=5; (${phr} / ${tot})*100" | bc -l)
echo -e "PolyHighResolution\t${phr}\t${perPHR}" >>  ${OUTDIR}/summary_${REF_ID}.txt
perNMH=$(echo "scale=5; (${nmh} / ${tot})*100" | bc -l)
echo -e "NoMinorHomt\t${mnh}\t${perMNH}" >>  ${OUTDIR}/summary_${REF_ID}.txt
perMHR=$(echo "scale=5; (${mhr} / ${tot})*100" | bc -l)
echo -e "MonoHighResolution\t${mhr}\tperMHR" >>  ${OUTDIR}/summary_${REF_ID}.txt
perOTHER=$(echo "scale=5; (${other} / ${tot})*100" | bc -l)
echo -e "Other\t${other}\t${perOTHER}" >>  ${OUTDIR}/summary_${REF_ID}.txt
perCRBT=$(echo "scale=5; (${crbt} / ${tot})*100" | bc -l)
echo -e "CallRateBelowThreshold\t${crbt}\t${perCRBT}" >>  ${OUTDIR}/summary_${REF_ID}.txt
perHEM=$(echo "scale=5; (${hem} / ${tot})*100" | bc -l)
echo -e "Hemizygous\t${hem}\t${perHEM}" >>  ${OUTDIR}/summary_${REF_ID}.txt
perOTV=$(echo "scale=5; (${otv} / ${tot})*100" | bc -l)
echo -e "OTV\t${otv}\t${perOTV}" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Number of SNPs: ${tot}" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo ""  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "Sample QC Thresholds" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "axiom_dishqc_DQC: >=  $DQC_THRESHOLD" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "qc_call_rate: >=  $CALL_RATE_THRESHOLD" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "plate_qc_percentSamplePassed: >=  $PLATE_PASS_RATE_THRESHOLD" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "plate_qc_averageCallRate: >=  $AVERAGE_PLATE_CALL_RATE" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo ""  >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "SNP QC Thresholds" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "##########################" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "specie-type: >=  ${SPECIE}" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "cr-cutoff: >=  $CR_CUTOFF" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "fld-cutoff: >=  $FLD_CUTOFF" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "het-so-cutoff: >=  $HET_SO_CUTOFF" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "het-so-otv-cutoff: >=  $HET_SO_OTV_CUTOFF" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "hom-ro-1-cutoff: >=  $HOM_RO_1_CUTOFF" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "hom-ro-2-cutoff: >=  $HOM_RO_2_CUTOFF" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "hom-ro-3-cutoff: >=  $HOM_RO_3_CUTOFF" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "hom-ro: $HOM_RO_BOOL" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "hom-het: $HOM_HET_BOOL" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "num-minor-allele-cutoff: >=  $NMA_CUTOFF" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "priority-order: $PRIORITY" >>  ${OUTDIR}/summary_${REF_ID}.txt
echo "recommended: $RECOMMENDED" >>  ${OUTDIR}/summary_${REF_ID}.txt



