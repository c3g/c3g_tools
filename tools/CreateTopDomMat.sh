#!/bin/env bash

if [ $# -ne 2 ]; then
    echo $0: usage: CreateMatA.sh myfile myRes
    exit 1
fi

input=$1
res=$2

## create file format
tail -n +2  $input > ${input}.tmp
awk -v res=${res} '{n=split($1, a, "-"); sub(/-[0-9]+$/, "", $1);print $1"\t"a[n]"\t"(a[n]+res)}' ${input}.tmp > ${input}.bins
cut -f 2- ${input}.tmp > ${input}.tmp2
paste ${input}.bins ${input}.tmp2 > ${input}.MatA
rm -f ${input}.tmp ${input}.tmp2 ${input}.bins
