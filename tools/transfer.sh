#!/bin/bash
set -eu -o pipefail

usage () {
echo "$0 <n_thread> <source> <dest>"
echo "a threaded copy tool" 
}

if [[ $# -ne 3 ]] ; then
 usage
fi

export SOURCE_DIR=$2
export DEST_DIR=$3
export threads=$1

mkdir -p  $DEST_DIR 2>/dev/null
# sync folder structure first
rsync -a -f'+ */' -f'- *' $SOURCE_DIR $DEST_DIR

# cwd
cd $SOURCE_DIR

find . -type f -print0 |  parallel -0 --linebuffer --jobs $threads 'rsync -av {} $DEST_DIR/{//}/'
