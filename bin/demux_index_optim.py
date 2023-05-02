#!/usr/bin/env python

import sys
import os
import argparse
import re
import gzip

#barcode_path='/Users/martingordon/Documents/projects/041223_Ronald_B_Babu/041323_RBenjamin_gigaassay/clustering_test/starcode_umi_clusters_index_10k.txt'
#fastq_path='/Users/martingordon/Documents/projects/041223_Ronald_B_Babu/041323_RBenjamin_gigaassay/clustering_test/SRR20707787_umi_10000.fastq.gz'
#sample_name='test-index'
#output_dir='./demux_index'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process fastq files and split reads into new files based on the barcode found at the 3\' end of the reads.')
    parser.add_argument('--fastq_path', metavar='fastq_path', type=str, help='path to the input fastq file')
    parser.add_argument('--barcode_path', metavar='barcode_path', type=str, help='path to the input barcode file')
    parser.add_argument('--sample_name', metavar='sample_name', type=str, help='output file name' )
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
    

#read in starcode umi clusters & represent as dictionary
#update: index as key & barcodes as val for faster searching
def read_barcodes(barcode_path):
    barcodes = {}
    with gzip.open(barcode_path, 'rt') as barcode_file:
        for line in barcode_file:
            barcode_group, barcode_count, barcode_index = line.strip().split('\t')
            for barcode in barcode_index.split(','): #split each index
                barcode = int(barcode) #make int for faster searching
                if barcode not in barcodes: #check each n is unique
                    barcodes[barcode] = set()
                barcodes[barcode].add(barcode_group)
    return barcodes

#match fq record index with barcode_dict key; return barcode consensus seq 
def match_index(fq_record, barcodes):
    barcode_group = None
    if fq_record in barcodes:
        barcode_group = barcodes[fq_record]
        del barcodes[fq_record]  # discard index to decrease future search space
    return str(barcode_group).split("'")[1]
    #return str(barcode_group)    


def process_fastq(fastq_path, barcode_path, output_dir, sample_name):
    barcodes = read_barcodes(barcode_path)

    with gzip.open(fastq_path,'rt',encoding='utf-8') as fastq_file:
        for count, (read_id, seq, plus, qual) in enumerate(zip(*[fastq_file]*4), 1): # zip to group the four lines of record, 1 indexing to match starcode
                barcode_group = match_index(count, barcodes)
                if barcode_group is not None and len(barcode_group) == 32: #only want reads with 32bp barcodes
                    output_file_name = os.path.join(output_dir, '%s_%s.demux.fastq.gz' % (barcode_group, sample_name))
                    with gzip.open(output_file_name, 'at', encoding='utf-8', compresslevel=6) as output_file:
                        output_file.write('%s\n%s\n%s\n%s\n' % (read_id.strip(), seq.strip(), plus.strip(), qual.strip()))
                else:
                    barcode_group = 'unassigned'
                    output_file_name = os.path.join(output_dir, '%s_%s.demux.fastq.gz' % (barcode_group, sample_name))
                    with gzip.open(output_file_name, 'at', encoding='utf-8', compresslevel=6) as output_file:
                        output_file.write('%s\n%s\n%s\n%s\n' % (read_id.strip(), seq.strip(), plus.strip(), qual.strip()))

#run
process_fastq(args.fastq_path,args.barcode_path,args.output_dir,args.sample_name)