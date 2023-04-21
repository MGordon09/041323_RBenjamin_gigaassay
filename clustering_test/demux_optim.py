#!/usr/bin/env python

import sys
import os
import argparse
import re
import gzip

def read_barcodes(barcode_path):
    barcodes = {}
    with gzip.open(barcode_path, 'rt') as barcode_file:
        for line in barcode_file:
            barcode_group, barcode_count, barcode_list = line.strip().split('\t')
            barcodes[barcode_group] = set(barcode_list.split(','))
    return barcodes

def match_barcode(read_seq, barcodes):
    for barcode_group, barcode_list in barcodes.items():
        for barcode in barcode_list:
            if read_seq[-32:].endswith(barcode):
                return barcode_group
    return "unknown"


def process_fastq(fastq_path, barcode_path, output_dir, sample_name):
    barcodes = read_barcodes(barcode_path)

    with gzip.open(fastq_path,'rt',encoding='utf-8') as fastq_file, gzip.open(output_dir + sample_name +'_unassigned.fastq', 'wt', encoding='utf-8') as unknown_file:
        for line in fastq_file:
            if line.startswith('@'):  # Header line for a read
                read_id = line.strip()
                seq = fastq_file.readline().strip()
                plus = fastq_file.readline().strip() 
                qual = fastq_file.readline().strip()

                barcode_group = match_barcode(seq, barcodes)
                
                output_file_name = os.path.join(output_dir, '%s_%s.demux.fastq' % (barcode_group, sample_name))
                with open(output_file_name, 'a') as output_file:
                    output_file.write('%s\n%s\n%s\n%s\n' % (read_id, seq, plus, qual))
                    
                if barcode_group == "unknown":
                    unknown_file.write('%s\n%s\n%s\n%s\n' % (read_id, seq, plus, qual))
                    
    return


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