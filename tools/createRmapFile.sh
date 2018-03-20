##!/bin/env sh
# 
##2018-03-15 - Edouard Henrion - edouard.henrion@computationalgenomics.ca
#Usage createRmapFile.sh genome_digest_file($1) output_rmap_file($2)

GENOME_DIGEST=$1
RMAP_FILE=$2

cut -f 1-3 $GENOME_DIGEST > ${RMAP_FILE}.tmp && \
awk 'BEGIN {{FS=\"\\t\"; OFS=\"\\t\"}} NR>2 {{if ($2 != $3) print $0, NR}}' ${RMAP_FILE}.tmp > $RMAP_FILE && \
rm ${RMAP_FILE}.tmp
