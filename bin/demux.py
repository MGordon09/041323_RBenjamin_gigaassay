#!/usr/bin/env python

import sys
import os
import argparse
import re


#TODO maybe trim final 32bp of reads from file or just do prior? maybe use cutadapt as optimised..
#take gzipped files as input and process for increased speed

def read_barcodes(barcode_path):
    barcodes = {}
    with open(barcode_path, 'r') as barcode_file:
        for line in barcode_file:
            barcode_group, barcode_count, barcode_list = line.strip().split('\t')
            barcodes[barcode_group] = set(barcode_list.split(','))
    return barcodes

def match_barcode(read_seq, barcodes):
    for barcode_group, barcode_list in barcodes.items():
        for barcode in barcode_list:
            if read_seq[-32:].endswith(barcode):
                return barcode_group #, read_seq[:-len(barcode))]? trim final 32bp from reads?
    return "unknown"

# def match_barcode(read_seq, barcodes):
#     for barcode_group, barcode_list in barcodes.items():
#         for barcode in barcode_list:
#             pattern = f"{barcode}$"
#             if re.search(pattern, read_seq[-32:]): #only search final 32 characters of read for match
#                 return barcode_group
#     return "unknown"


def process_fastq(fastq_path, barcode_path, output_dir, sample_name):
    barcodes = read_barcodes(barcode_path)
    barcode_counts = {}

    with open(fastq_path, 'r') as fastq_file, open(output_dir + sample_name +'_unassigned.fastq', 'w') as unknown_file:
        for line in fastq_file:
            if line.startswith('@'):  # Header line for a read
                read_id = line.strip()
                seq = fastq_file.readline().strip()
                plus = fastq_file.readline().strip() 
                qual = fastq_file.readline().strip()

                barcode_group = match_barcode(seq, barcodes)
                barcode_counts[barcode_group] = barcode_counts.get(barcode_group, 0) + 1 #counter for reads beloning to barcode group
                
                output_file_name = os.path.join(output_dir, '%s_%s.demux.fastq' % (barcode_group, sample_name))
                with open(output_file_name, 'a') as output_file:
                    output_file.write('%s\n%s\n%s\n%s\n' % (read_id, seq, plus, qual))
                    
                if barcode_group == "unknown":
                    unknown_file.write('%s\n%s\n%s\n%s\n' % (read_id, seq, plus, qual))
                    
    return barcode_counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process fastq files and split reads into new files based on the barcode found at the 3\' end of the reads.')
    parser.add_argument('--fastq_path', metavar='fastq_path', type=str, help='path to the input fastq file')
    parser.add_argument('--barcode_path', metavar='barcode_path', type=str, help='path to the input barcode file')
    parser.add_argument('--sample_name', metavar='sample_name', type=str, help='output file name')
    parser.add_argument('--output_dir', metavar='output_dir', type=str, help='output directory (default: current directory)', default='.')
    args = parser.parse_args()
    
    if not os.path.exists(args.fastq_path):
        raise ValueError(f'File not found: {args.fastq_path}')
    if not os.path.exists(args.barcode_path):
        raise ValueError(f'File not found: {args.barcode_path}')
    if args.sample_name is None:
        raise ValueError(f'Output name not given: {args.sample_name}')
    if not os.path.isdir(args.output_dir):
        raise ValueError(f'Output directory not found: {args.output_dir}')

    barcode_counts = process_fastq(args.fastq_path, args.barcode_path, args.output_dir, args.sample_name)
    
    print(f'Total reads processed: {sum(barcode_counts.values())}')
    print(f'Total barcode groups found: {len(barcode_counts)}')
    print('Barcode group counts:')
    for barcode_group, count in barcode_counts.items():
        print(f'{barcode_group}: {count}')