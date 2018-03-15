#!/bin/env bash
# Exit immediately on error
#set -eu -o pipefail

################################################################################
# This program prepares the data needed to run pipeline tests and show how to
# run the standard pipeline for the first time using the default parameters
################################################################################

create_dir() {
  DIR=$1
  # Create directory with permissions if necessary
  if [[ ! -d $DIR ]]
  then
    mkdir -p $DIR
    chmod ug+rwX,o+rX-w $DIR
  fi
}

check_files(){
  # if design file is in the DATA_PATH directory, take this as the complete path to design file
  if [[ -f  $DATA_PATH/$DESIGN_FILE ]]; then
    DESIGN_PATH=$DATA_PATH/$DESIGN_FILE
  else
    DESIGN_PATH=$DESIGN_FILE
  fi  
  # check the readset File 
  if [[ -f  $DATA_PATH/$READSET_FILE ]]; then
    READSET_PATH=$DATA_PATH/$READSET_FILE
  else
    READSET_PATH=$DATA_PATH/$READSET_FILE
  fi
  # check the target bed 
  if [[ -f  $DATA_PATH/$BED_FILE ]]; then
    BED_PATH=$DATA_PATH/$BED_FILE
  else
    BED_PATH=$DATA_PATH/$BED_FILE
  fi
  if [[ ! -f  $DESIGN_PATH ]]; then
     echo "" >&2
     echo "Design file not found or not specified: "$DESIGN_PATH", defaults to design.tsv in:" $OUTPUT_DIR >&2
     DESIGN_PATH=$OUTPUT_DIR/design.tsv 
     echo "" >&2
  else
     cp $DESIGN_PATH $OUTPUT_DIR/
  fi
  if [[ ! -f  $READSET_PATH ]]; then
     echo "" >&2
     echo "Readset file not found or not specified: "$READSET_PATH", defaults to readset.tsv in:" $OUTPUT_DIR  >&2
     echo "" >&2
     READSET_PATH=$OUTPUT_DIR/readset.tsv 
  else
    cp $READSET_PATH $OUTPUT_DIR/
  fi
  if [[ ! -f  $BED_PATH ]]; then
     echo "" >&2
     echo "BED file not found or not specified: "$BED_PATH" in:" $OUTPUT_DIR >&2
     echo "" >&2
     BED_PATH=$OUTPUT_DIR/ 
  else
    cp $BED_PATH $OUTPUT_DIR
  fi
}


# synchronize files from repository
copy_files(){
    # Get sample names, create the include command for rsync
    INCLUDE_RAW=`cat $READSET_PATH | sed '1d' | awk -F"\t" '{print $1}' | uniq -u | awk -F"\t" '{print "--include \x27raw_reads/"$1"/\x27 --include \x27raw_reads/"$1"/*\x27"}'`
    INCLUDE_TRIM=`cat $READSET_PATH | sed '1d' | awk -F"\t" '{print $1}' | uniq -u | awk -F"\t" '{print "--include \x27trim/"$1"/\x27 --include \x27trim/"$1"/*\x27"}'`
    # rsync symbolic links as symbolic links
    RSYNC_COMMAND=`echo "rsync -avPl --include raw_reads/ --include trim/ "$INCLUDE_RAW" "$INCLUDE_TRIM" --exclude='*' "$DATA_PATH"/ "$OUTPUT_DIR"/"`   
    echo $RSYNC_COMMAND >&2
    eval $RSYNC_COMMAND
}


# Install the latest version of mugqic_pipelines in $INSTALL_DOWNLOAD 
get_mugqic_pipelines_zip(){
  SOFTWARE=mugqic_pipelines
  ARCHIVE=$SOFTWARE-$VERSION.tar.gz
  ARCHIVE_URL=https://bitbucket.org/mugqic/$SOFTWARE/downloads/$ARCHIVE  
  curDir=`pwd`
  cd $INSTALL_DOWNLOAD
  wget --no-check-certificate $ARCHIVE_URL --output-document=$INSTALL_DOWNLOAD/$ARCHIVE
  tar zxvf $ARCHIVE
  echo "" >&2
  echo "Mugqic pipelines installed in:" $INSTALL_DOWNLOAD >&2
  echo "" >&2
  cd $curDir
  PIPELINE_DIR=`echo $INSTALL_DOWNLOAD/$SOFTWARE-$VERSION/pipelines/$PIPELINE`
}

# If a specific commit is required, the installation will change slightly
get_mugqic_pipelines_commit(){
  SOFTWARE=mugqic_pipelines
  VERSION=$1
  ARCHIVE=$SOFTWARE-$VERSION.tar.gz
  ARCHIVE_URL=https://bitbucket.org/mugqic/$SOFTWARE/downloads/$ARCHIVE  
  curDir=`pwd`
  cd $INSTALL_DOWNLOAD
  if [[ -e $INSTALL_DOWNLOAD/mugqic_pipelines ]]; then
    echo "$INSTALL_DOWNLOAD/mugqic_pipelines already exist; it will be updated!" >&2
  else
    git clone git@bitbucket.org:mugqic/mugqic_pipelines.git
  fi
  cd mugqic_pipelines
  git checkout $VERSION ./
  echo "" >&2
  echo "Mugqic pipelines commit "$VERSION" installed in:" $INSTALL_DOWNLOAD >&2
  echo "" >&2
  cd $curDir
  PIPELINE_DIR=`echo $INSTALL_DOWNLOAD/mugqic_pipelines/pipelines/$PIPELINE` 
}

# Copy from a local directory
get_mugqic_pipelines_local(){	
    VERSION=$1
    curDir=`pwd`
    cd $INSTALL_DOWNLOAD
    if [[ -e $INSTALL_DOWNLOAD/mugqic_pipelines ]]; then
        echo "$INSTALL_DOWNLOAD/mugqic_pipelines already exist; it will be updated!" >&2
    else
        mkdir -p $INSTALL_DOWNLOAD/mugqic_pipelines
    fi	
    RSYNC_COMMAND=`echo "rsync -avPL "$VERSION"/mugqic_pipelines "$INSTALL_DOWNLOAD""`   
    echo $RSYNC_COMMAND >&2
    eval $RSYNC_COMMAND     
    echo "" >&2
    echo "Mugqic pipelines copied from local directory "$VERSION" installed in:" $INSTALL_DOWNLOAD >&2
    echo "" >&2
    cd $curDir
    PIPELINE_DIR=`echo $INSTALL_DOWNLOAD/mugqic_pipelines/pipelines/$PIPELINE` 
}

# Call pipeline
call_pipeline(){
  echo $PIPELINE_DIR/$PIPELINE.py "-h" 
  TOSTEP=
  #TOSTEP=`eval $PIPELINE_DIR/$PIPELINE.py  2>&1 | grep '^[0-9]*[\-].*$' | awk -F"-" 'BEGIN {max = 0} {if ($1>max) max=$1} END {print max}'`
  if [[ -f  $DESIGN_PATH ]]; then 
    echo -e "$PIPELINE_DIR/$PIPELINE.py -c $PIPELINE_DIR/$PIPELINE.base.ini -s 1-$TOSTEP -d $DESIGN_FILE -r $READSET_FILE > toRun.sh"
  else
    echo -e "$PIPELINE_DIR/$PIPELINE.py -c $PIPELINE_DIR/$PIPELINE.base.ini -s 1-$TOSTEP -r $READSET_FILE > toRun.sh"
  fi
}

DESIGN_PATH=
READSET_PATH=
BED_PATH=
PIPELINE_CALL=
PIPELINE_DIR=

create_dir $OUTPUT_DIR
cd $OUTPUT_DIR
create_dir $INSTALL_DOWNLOAD

if [[ "$VERSION" =~ [^0-9]*[\.].* ]]; then
  echo "NOTICE: will install pipeline using a version number: "$VERSION >&2
  get_mugqic_pipelines_zip
elif [[ -d "$VERSION" ]]; then
    echo "NOTICE: will install pipeline from a local copy : "$VERSION >&2
    get_mugqic_pipelines_local $VERSION
else
  echo "NOTICE: will install pipeline using a commit number: "$VERSION >&2
  get_mugqic_pipelines_commit $VERSION
fi

# Check design, readset and bed files
check_files
# Copy files using rsync 
echo "NOTICE: copying (rsync) files from "$DATA_PATH " to "  $OUTPUT_DIR >&2
copy_files
# Call the pipeline and generates the first command
echo "## NOTICE: Pipeline usage" 
call_pipeline