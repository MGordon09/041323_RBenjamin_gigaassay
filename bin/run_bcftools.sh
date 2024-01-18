#!/usr/bin/env bash

# ------------------------------------------------------------------
#  Run bcftools mpileup and variant calling on each file 
#
#set -e  # Abort script at first error, when a command exits with non-zero status (except in until or while loops, if-tests, list constructs)
#set -u  # Attempt to use undefined variable outputs error message, and forces an exit
#set -x  # Similar to verbose mode (-v), but expands commands
#set -o pipefail  # Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.

# parameters
REFERENCE=$1
BAM=$2

# extract sample id and barcode from read name
FILE="$(basename $BAM)"

SAMPLE_ID="${FILE%%_*}" #extract the barcode from file name
BARCODE_ID="${FILE#*_}" # extract everything after _ and strip extension
BARCODE_ID="${BARCODE_ID%%.*}"

bcftools mpileup \
    -I -Ou \
    -f $REFERENCE \
    $BAM \
| bcftools call \
    --ploidy 1 \
    --skip-variants indels \
    --consensus-caller \
    --variants-only \
    -Ou -o ./variant.files/${SAMPLE_ID}_${BARCODE_ID}.bcf 2> ./variant.files/${SAMPLE_ID}_${BARCODE_ID}.bcf.err