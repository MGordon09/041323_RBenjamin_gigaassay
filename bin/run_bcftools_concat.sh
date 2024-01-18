#!/usr/bin/env bash

# ------------------------------------------------------------------
#  Run bcftools concat on the batch list of files
# 
#
#set -e  # Abort script at first error, when a command exits with non-zero status (except in until or while loops, if-tests, list constructs)
#set -u  # Attempt to use undefined variable outputs error message, and forces an exit
#set -x  # Similar to verbose mode (-v), but expands commands
#set -o pipefail  # Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.

# parameters
SAMPLE_ID=$1
VCF_LIST=$2

# extract the vcf subset name
BATCH="$(basename $VCF_LIST)"

echo "combining $BATCH, sorting..."

bcftools concat \
    --allow-overlaps \
    --file-list ${VCF_LIST} \
    | bcftools sort -o ${SAMPLE_ID}.${BATCH}.subset.vcf


bcftools view -Oz -o ${SAMPLE_ID}.${BATCH}.subset.vcf.gz ${SAMPLE_ID}.${BATCH}.subset.vcf

echo "Indexing ${BATCH}..."

# index intermediate files to merge
bcftools index ${SAMPLE_ID}.${BATCH}.subset.vcf.gz