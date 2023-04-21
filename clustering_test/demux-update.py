import sys
import os
import gzip

"""
Read a fastq (or fastq.gz) file, split according to the barcodes given in
files in the format output by starcode software (https://github.com/gui11aume/starcode)
"""

fq_f = sys.argv[1]
bc_f = sys.argv[2]

#create dict from barcode file
barcode_dict = {}

with open(bc_f) as barcode_file:
    for line in barcode_file:
        barcode, _, seqs = line.rstrip().split('\t')
        barcode_dict[barcode] = seqs.split(',')

#create output files for each barcode
out_files = {}
for bar in barcode_dict.keys():
    fname = f'./{bar}_placeholder_.fastq'
    out_files[bar] = open(fname, 'w')

#output file for unknown barcodes
unknown_fname = f'./unknown_placeholder_.fastq'
unknown_file = open(unknown_fname, 'w')

#read the fastq file
with gzip.open(fq_f, 'rt') if fq_f.endswith('.gz') else open(fq_f) as fastq_file:
    while True:
        header = fastq_file.readline().rstrip()
        if not header:
            break # end of file
        sequence = fastq_file.readline().rstrip()
        plus = fastq_file.readline().rstrip()
        quality = fastq_file.readline().rstrip()

        #find the matching barcode
        detected = False
        for barcode, bar_seqs in barcode_dict.items():
            for b in bar_seqs:
                if sequence.endswith(b):
                    detected = True
                    out_files[barcode].write(f'{header}\n{sequence}\n{plus}\n{quality}\n')
                    break # exit the loop if match found
            if detected:
                break # exit the loop if match found

        #write to the unknown file if no match found
        if not detected:
            unknown_file.write(f'{header}\n{sequence}\n{plus}\n{quality}\n')

#close all the output files
for out_file in out_files.values():
    out_file.close()
unknown_file.close()