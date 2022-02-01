#!/bin/bash
set -eu -o pipefail

usage () {
echo "$0 <n_thread> <source> <dest>"
echo "a threaded copy tool that will rsync <source> inside <dest>"
}

if [[ $# -ne 3 ]] ; then
 usage
fi

export SOURCE_DIR=$2
export DEST_DIR=$3
export threads=$1

export DEST_DIR=${DEST_DIR%/}

# mkdir -p  $DEST_DIR 2>/dev/null
# sync folder structure first
rsync -L -a -f'+ */' -f'- *' $SOURCE_DIR $DEST_DIR

# cwd
cd $SOURCE_DIR/..
SOURCE_DIR_LAST=`basename $SOURCE_DIR`
find $SOURCE_DIR_LAST/ -type f -print0 |  parallel -0 --linebuffer --jobs $threads 'rsync -L -av {} $DEST_DIR/{//}/'
