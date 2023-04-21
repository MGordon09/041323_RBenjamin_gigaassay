import sys
import os
import gzip

"""
Read a fastq (or fastq.gz) file, split according to the barcodes given in
files in the format output by starcode software (https://github.com/gui11aume/starcode)
"""

fq_f=sys.argv[1]
bc_f=sys.argv[2]


#create dict from barcode file
bases=("A","G","C","T","U")
barcode_dict = {}

with open(bc_f) as barcode_file:
    for line in barcode_file:
        cols=line.split('\t')
        barcode=cols[0]
        seqs=cols[2]
        barcode_dict[barcode] = seqs

#iterate through the keys
for bar in barcode_dict.keys():

    #list of barcodes w errors
    bar_seqs=barcode_dict.get(bar).split(',')

    #output filename
    fname='./' + str(bar) + '_' + 'placeholder' + '_' '.fastq' 
    
    with open(fq_f) as fastq_file, open(fname, 'w') as out_file:
        while True:
            header=fastq_file.readline().strip()
            if not header:
                break # end of file
            sequence=fastq_file.readline().strip()
            plus=fastq_file.readline().strip()
            quality=fastq_file.readline().strip()

            # if match found:
            detected=False

            for b in bar_seqs:
                if sequence.endswith(b):
                    print(b)
                
                    detected=True

                    out_file.write(header + '\n')
                    out_file.write(sequence + '\n')
                    out_file.write(plus + '\n')
                    out_file.write(quality + '\n')

            if not detected:
                unknown='./' + 'unknown' + '_' + 'placeholder' + '_' '.fastq'
                with open(unknown, 'w') as unknown_file:
                    unknown_file.write(header + '\n')
                    unknown_file.write(sequence + '\n')
                    unknown_file.write(plus + '\n')
                    unknown_file.write(quality + '\n')