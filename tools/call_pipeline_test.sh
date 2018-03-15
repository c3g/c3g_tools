#!/bin/env sh
# Exit immediately on error
#set -eu -o pipefail

################################################################################
# This is a pipeline test program that should be used to shortcut the 
# use of the mugqic pipeline tests templates. Use this program
# to ensure the consistency between pipeline tests
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

usage() {
SCRIPT=`basename ${BASH_SOURCE[0]}`
usage="\
Usage: $SCRIPT [OPTION]... -r [DATA_PATH] -R [READSET_FILE] -d [DESIGN_FILE] -b [BED_FILE] -o [OUTPUT_DIR] -i [INSTALL_DOWNLOAD] -p [PIPELINE] -V [VERSION]

Options:
-h  display this help and exit.
-r    DATA_PATH   Path to the data repository 
-R    READSET_FILE   Name of the design file (could be set as a complete path or a file name, if file name is set, the file must be located in DATA_PATH)
-d    DESIGN_FILE   Name of the design file (could be set as a complete path or a file name, if file name is set, the file must be located in DATA_PATH)
-b    BED_FILE    Name of the readsets target coverage BED file (could be set as a complete path or a file name, if file name is set, the file must be located in DATA_PATH)
-o    OUTPUT_DIR  Destination folder, default is current directory
-i    INSTALL_DOWNLOAD    Pipeline install download PATH, default is OUTPUT_DIR/bin
-p    PIPELINE    pipeline to run   (one of the following: dnaseq | rnaseq | chipseq | rnaseq_denovo_assembly )
-V    VERSION   pipeline version or commit number (by default the last commit in the master repository is installed)

Environment variables override the default commands:
DATA_PATH DESIGN_FILE READSET_FILE BED_FILE OUTPUT_DIR INSTALL_DOWNLOAD PIPELINE VERSION
"
    echo $usage
} 

# Destination folder, default is current directory
OUTPUT_DIR=`pwd`
## Install download, default is $OUTPUT_DIR/bin
INSTALL_DOWNLOAD=$OUTPUT_DIR/bin 
DESIGN_FILE=
READSET_FILE=
BED_FILE=
VERSION="master"
DATA_PATH=
PIPELINE=


NUMARGS=$#
echo -e \\n"Number of arguments: $NUMARGS"
if [ $NUMARGS -eq 0 ]; then
  usage
fi

while getopts "p:V:r:R:d:b:o:i:h" OPTION; 
do
  case $OPTION in
    h) 
    usage
    ;;
    r)
    DATA_PATH=$OPTARG 
    ;;
    d)
    DESIGN_FILE=$OPTARG 
    ;;
    R)
    READSET_FILE=$OPTARG 
    ;;
    b)
    BED_FILE=$OPTARG 
    ;;
    o)
    OUTPUT_DIR=$OPTARG 
    ;;
    i)
    INSTALL_DOWNLOAD=$OPTARG
		if  [[ -z $INSTALL_DOWNLOAD ]] ; then 
			INSTALL_DOWNLOAD=$OUTPUT_DIR/bin 
		fi
		;;
    p)
    PIPELINE=$OPTARG
        case $PIPELINE in
          dnaseq | rnaseq | chipseq | rnaseq_denovo_assembly ) echo "$0: Pipeline to run : $PIPELINE" >&2 ;;
          *) echo "$0: Invalid pipeline : $PIPELINE" >&2
            ;; 
        esac  
    ;;
   V)
   VERSION=$OPTARG
	;;
   ?)
   usage ;;
  esac
done


echo -e "calling pipeline tests with parameters -r $DATA_PATH -p $PIPELINE -V $VERSION -R $READSET_FILE -d $DESIGN_FILE -o $OUTPUT_DIR -i $INSTALL_DOWNLOAD -b $BED_FILE"
# Stop if mandatory variables are not set
if [[ -n $DATA_PATH ]] &&  [[ -n $PIPELINE ]]  && [[  -n $VERSION ]] ; then
  # Call generic module install script once all variables and functions have been setup
  PIPELINE_TEST_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  source $PIPELINE_TEST_SCRIPT_DIR/pipeline_test.sh $@
else
 if  [[ -z $DATA_PATH ]] ; then 
    echo "$0: Path to the data repository is mandatory: $DATA_PATH" >&2
 elif  [[ -z $PIPELINE ]] ; then 
    echo "$0: Pipeline is mandatory: $PIPELINE" >&2
 elif  [[ -z $VERSION ]]; then 
     echo "$0: Pipeline Version is mandatory: $VERSION" >&2
 else
    echo "$0: Use the mandatory parameters: $DATA_PATH $PIPELINE $VERSION" >&2
 fi
fi