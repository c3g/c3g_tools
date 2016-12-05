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

exit 0

}


##Default args
APT_MODULE=mugqic_dev/AffymetricxApt/1.18.0
MUGIC_TOOLS_MODULE=mugqic_dev/mugqic_tools/2.1.7
#r module should have the SNPolisher package installed
R_MODULE=mugqic_dev/R_Bioconductor/3.2.3_3.2
DQC_THRESHOLD=0.82
CALL_RATE_THRESHOLD=97
PLATE_PASS_RATE_THRESHOLD=95
AVERAGE_PLATE_CALL_RATE=98.5
SPECIE=human
OUTPUT_SNP_NUMBER=6

if [ $# -lt 12 ]
then
	usage
	exit 0
fi

##get arguments
while getopts "A:M:R:c:l:o:m:a:r:d:s:p:C:S:n:" OPT
do
    case "$OPT" in
        A) 
           APT_MODULE=$OPTARG
           ;;
	M) 
           MUGIC_TOOLS_MODULE=$OPTARG
           ;;
        R)
           R_MODULE=$OPTARG
           ;;
        c)
           CEL_PATH=$OPTARG
           ;;
        l)
           ANALYSIS_FILES_DIR=$OPTARG
           ;;
        o)
           OUTDIR=$OPTARG
           ;;
        m)
           MASTER_LIST=$OPTARG
           ;;
        a)
           AXIOM_ARRAY_NAME=$OPTARG
           ;;
        r)
           AXIOM_ARRAY_REV=$OPTARG
           ;;
        d)
           DQC_THRESHOLD=$OPTARG
           ;;
        s)
          CALL_RATE_THRESHOLD=$OPTARG
           ;;
        p)
           PLATE_PASS_RATE_THRESHOLD=$OPTARG
           ;;
        C)
           AVERAGE_PLATE_CALL_RATE=$OPTARG
           ;;
        
        S) 
           SPECIE=$OPTARG
           ;;
        n)
           OUTPUT_SNP_NUMBER=$OPTARG
           ;;
        
    esac

done


REF_ID=$(date +"%F_%H-%M-%S")

#load APT 
module load ${APT_MODULE} ${MUGIC_TOOLS_MODULE} ${R_MODULE}

## dev_argument
#R_TOOLS=/lb/project/mugqic/projects/GeneTitan_BFXTD30/scripts/


mkdir -p ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/  ${OUTDIR}/GENO_QC_${REF_ID}/

########
## step1 Group Samples into Batches
## create file containing the list of cel files path
########
echo "-----------------------"
echo "Runing Step 1 - Sample Grouping"
echo "-----------------------"

echo cel_files > ${OUTDIR}/cel_list1_${REF_ID}.txt
ls ${CEL_PATH}/*.CEL >>  ${OUTDIR}/cel_list1_${REF_ID}.txt

continue=$(wc -l ${OUTDIR}/cel_list1_${REF_ID}.txt | cut -d\  -f1)

if [[ $continue == 1 ]]
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

continue=$(wc -l ${OUTDIR}/cel_list2_${REF_ID}.txt | cut -d\  -f1)

if [[ $continue == 1 ]]
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

continue=$(wc -l ${OUTDIR}/cel_list3_${REF_ID}.txt | cut -d\  -f1)

if [[ $continue == 1 ]]
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
-c ${OUTDIR}/cel_list3_${REF_ID}.txt \
-d ${PLATE_PASS_RATE_THRESHOLD} \
-a ${AVERAGE_PLATE_CALL_RATE} \
-p ${MASTER_LIST} \
-o ${OUTDIR}/cel_list4_${REF_ID}.txt"

Rscript ${R_TOOLS}/filterAxiom.R -s plateqc \
-q ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/AxiomGT1.report.txt \
-c ${OUTDIR}/cel_list3_${REF_ID}.txt \
-d ${PLATE_PASS_RATE_THRESHOLD} \
-a ${AVERAGE_PLATE_CALL_RATE} \
-p ${MASTER_LIST} \
-o ${OUTDIR}/cel_list4_${REF_ID}.txt


continue=$(wc -l ${OUTDIR}/cel_list4_${REF_ID}.txt | cut -d\  -f1)

if [[ $continue == 1 ]]
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
--ps2snp-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.ps2snp_map.ps"


ps-classification --species-type ${SPECIE} \
--metrics-file ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/PS_metrics.txt \
--output-dir ${OUTDIR}/GENOTYPE_AXIOM_${REF_ID}/ \
--ps2snp-file ${ANALYSIS_FILES_DIR}/${AXIOM_ARRAY_NAME}.${AXIOM_ARRAY_REV}.ps2snp_map.ps

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



