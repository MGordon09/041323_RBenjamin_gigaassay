#!/usr/bin/env bash

# ------------------------------------------------------------------
#  Run bbtools program demuxbyname.sh to 
#
set -e  # Abort script at first error, when a command exits with non-zero status (except in until or while loops, if-tests, list constructs)
set -u  # Attempt to use undefined variable outputs error message, and forces an exit
set -x  # Similar to verbose mode (-v), but expands commands
set -o pipefail  # Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.

# parameters
RESOURCES=$1
SAMPLES=$2
READS=$3

# extract sample id and barcode from read name
FILE="$(basename $READS)"
BARCODE_ID="${FILE%%_*}" #extract the barcode from file name

demuxbyname.sh \
    -Xmx=$RESOURCES \
    in=$READS \
    suffixmode=t \
    out=%_demux/${BARCODE_ID}.fastq.gz \
    names=$SAMPLES

# mv the data to subfolder
#while read s; do
#  mv ${BARCODE_ID}_${s}_demuxsamples.fastq.gz ./${s}_demux
#done < $SAMPLES