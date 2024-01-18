#!/usr/bin/env bash

# ------------------------------------------------------------------
#  Run bwa mem align on each PE file in folder
#
#set -e  # Abort script at first error, when a command exits with non-zero status (except in until or while loops, if-tests, list constructs)
#set -u  # Attempt to use undefined variable outputs error message, and forces an exit
#set -x  # Similar to verbose mode (-v), but expands commands
#set -o pipefail  # Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.

# parameters
RESOURCES=$1
READS=$2

# extract sample id and barcode from read name
FILE="$(basename $READS)"
SAMPLE_ID="${FILE%%_*}" #extract the sample name
BARCODE_ID="${FILE#*_}" # extract everything after _ and strip extension
BARCODE_ID="${BARCODE_ID%%.*}" # remove file suffix to 

# bwa mem index
INDEX_FLDR=$(find -L ./ -name "*.amb" | sed 's/.amb$//')

bwa mem \
    -t $RESOURCES -R "@RG\\tID:${BARCODE_ID}\\tSM:${SAMPLE_ID}\\tPL:Illumina" \
    $INDEX_FLDR 2> ./bam.files/${SAMPLE_ID}_${BARCODE_ID}.bwa.err $READS \
    | samtools sort -@ $RESOURCES -O bam -o ./bam.files/${SAMPLE_ID}_${BARCODE_ID}.sorted.bam