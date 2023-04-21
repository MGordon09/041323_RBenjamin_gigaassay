#!/usr/bin/env python

import sys
import os
import re
"""
Read a list of SE or PE fastq (or fastq.gz) files, split them according to the barcodes given in
files in the format
"""

## TODO 
#input files: use oslistdir to produce list of input files matching correct pattern
#include line to gzip input, or else just use process substitution




barcode_dict = {}

with open('./starcode_umi_toy.csv') as barcode_file:
    for line in barcode_file:
        cols=line.split('\t')
        barcode=cols[0]
        seqs=cols[2]
        barcode_dict[barcode] = seqs

#loop through input fastq files
# with open('./SRR20707787_trim_s1000.fastq.gz') as fastq_file:
#     for line in fastq_file:
#         if line


bases=("A","G","C","T","U")

# loop through consensus barcodes (KEYS)
for bar in barcode_dict.keys():
    # seqs is list of error-barcodes
    seqs=barcode_dict.get(bar).split(',') #list of barcodes
    #print('key: ' + str(bar))
    print('barcodes: '+ str(seqs))

    #if len(i) != 32: maybe think of filtering shorter barcodes?
    
    fname = "./" + str(bar) + '_' + 'samplename' + '_' '.fastq'
    #this will be the file to write output to 
    with open(fname, 'w') as file:    
        #loop through fastq files... 
        with open('./SRR20707787_trim_s1000.fastq') as fastq_file: #this will have to be list of files...
            for line in fastq_file:
                if line.startswith((bases)) and line.endswith(bar): #startswith can take tuple
                    print(line) 