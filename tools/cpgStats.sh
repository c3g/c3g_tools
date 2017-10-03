#!/usr/bin/env sh
set -e
set -o pipefail

# Get args
INPUT_CPG_PROFILE=$1
CG_STAT_OUTFILE=$2
LAMBDA_STAT_OUTFILE=$3
PUC19_STAT_OUTFILE=$4

### Calculate the stats for regular CpGs
sed 1d $INPUT_CPG_PROFILE | grep lambda -v | grep pUC19 -v | awk '{ cg1x+=1; if($11>=30) cg30x+=1; if($11>=10) cg10x+=1; }END {print cg1x","cg10x","cg30x}' > $CG_STAT_OUTFILE

# Calculate lambda conversion rate
cat $INPUT_CPG_PROFILE | grep lambda | awk '{meth+=$10;reads+=$11;} END {print (1-meth/reads)*100}' > $LAMBDA_STAT_OUTFILE

### Calculate the stats for pUC19
cat $INPUT_CPG_PROFILE | grep pUC19 | awk '$2<667 && $2>496 {meth+=$10;reads+=$11;} END {if (reads==0) print "NA"; else print (1-meth/reads)*100}' > $PUC19_STAT_OUTFILE

echo 0;
