#!/bin/env sh
set -e
set -o pipefail

OTU=$1
TMP="_temp"
BKP=$(echo $OTU | sed -e 's/\.txt$/_BACKUP.txt/')

grep -v "^[a-zA-Z]..*$" $OTU > $TMP && \
cp $OTU $BKP && \
mv $TMP $OTU
