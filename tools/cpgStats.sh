#!/usr/bin/env sh
set -o pipefail

# Get args
INPUT_CPG_PROFILE=$1
CG_STAT_OUTFILE=$2
LAMBDA_STAT_OUTFILE=$3
PUC19_STAT_OUTFILE=$4

echo -e "\n----------------------------------"
echo " CMD: cpgStats.sh  $INPUT_CPG_PROFILE $CG_STAT_OUTFILE $LAMBDA_STAT_OUTFILE $PUC19_STAT_OUTFILE"
echo -e "\nstep 1 regular CpGs "

### Calculate the stats for regular CpGs
sed 1d $INPUT_CPG_PROFILE | grep lambda -v | grep pUC19 -v | awk 'BEGIN{cg1x=0; cg10x=0; cg30x=0;}{ cg1x+=1; if($11>=30) cg30x+=1; if($11>=10) cg10x+=1; }END {print cg1x","cg10x","cg30x}' > $CG_STAT_OUTFILE
echo "...DONE"

# Calculate lambda conversion rate
echo  -e  "\nstep 2 Lambda conversion"
cat $INPUT_CPG_PROFILE | grep lambda  | awk ' BEGIN {meth=0; reads=0} {meth+=$10;reads+=$11} END {if (reads==0) {print "NA"} else {print (1-meth/reads)*100}}' > $LAMBDA_STAT_OUTFILE
echo "...DONE"


### Calculate the stats for pUC19
echo -e "\nstep 3 pUC19"
cat $INPUT_CPG_PROFILE | grep pUC19 | awk ' BEGIN {meth=0; reads=0} $2<667 && $2>496 {meth+=$10;reads+=$11} END {if (reads==0) {print "NA"} else {print (1-meth/reads)*100}}' > $PUC19_STAT_OUTFILE
echo "...DONE"

echo "----------------------------------"

exit
