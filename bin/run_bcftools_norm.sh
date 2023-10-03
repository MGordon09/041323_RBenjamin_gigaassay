#!/usr/bin/env bash

# ------------------------------------------------------------------
#  Run bcftools norm on the output files
# 
#
#set -e  # Abort script at first error, when a command exits with non-zero status (except in until or while loops, if-tests, list constructs)
#set -u  # Attempt to use undefined variable outputs error message, and forces an exit
#set -x  # Similar to verbose mode (-v), but expands commands
#set -o pipefail  # Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.

# parameters
REFERENCE=$1
VCF=$2

# extract sample id and barcode from read name
FILE="$(basename $VCF)"

SAMPLE_ID="${FILE%%_*}" #extract the barcode from file name
BARCODE_ID="${FILE#*_}" # extract everything after _ and strip extension
BARCODE_ID="${BARCODE_ID%%.*}"

if bcftools view  $VCF | grep -qv '^#'
then
    bcftools view -q 0.6:nref $VCF -Ou | \
    bcftools annotate --set-id "${SAMPLE_ID}_${BARCODE_ID}" |  \
    bcftools norm --fasta-ref $REFERENCE -Ob -o ${SAMPLE_ID}_${BARCODE_ID}.norm.bcf

    #create tbi index of file
    bcftools index ${SAMPLE_ID}_${BARCODE_ID}.norm.bcf
fi