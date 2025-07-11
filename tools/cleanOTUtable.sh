#!/bin/env sh
set -e
set -o pipefail

OTU=$1
TMP=$(mktemp)
BKP=$(echo $OTU | sed -e 's/\.txt$/_BACKUP.txt/')

grep -v "^[a-zA-Z]..*$" $OTU > $TMP && \
ls $(dirname $TMP) > /dev/null && \
cp $OTU $BKP && \
mv $TMP $OTU
