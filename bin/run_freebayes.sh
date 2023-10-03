#!/usr/bin/env bash

# ------------------------------------------------------------------
#  Run freebayes variant caller on each vcf file
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

freebayes \
    --min-alternate-count 1 \
    --min-alternate-fraction 0.6 \
    --min-mapping-quality 20  \
    --min-base-quality 20 \
    --use-best-n-alleles 3 \
    -f $REFERENCE \
    $BAM > ${SAMPLE_ID}_${BARCODE_ID}.vcf 2> ${SAMPLE_ID}_${BARCODE_ID}.vcf.err