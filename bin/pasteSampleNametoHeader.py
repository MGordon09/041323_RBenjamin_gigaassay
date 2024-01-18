#!/usr/bin/env python

# script to read in gzipped file, add sample name to end of read header and write out gzipped file

import sys
import os
import argparse
import gzip

# add arguments to the script
parser = argparse.ArgumentParser(description='Add sample names to read headers')

parser.add_argument(
    '--fastq_path', type = str, required = True,
    help = 'path to fastq file (gzipped)')

parser.add_argument(
    '--sample_name', type = str, required = True,
    help = "sample name/identifer to be added to read header")
    
parser.add_argument(
    '--output_dir', type = str, required = True,
    help = "output directory to write file to")

args = parser.parse_args()


#fastq_path = '/Users/martingordon/Documents/projects/041223_Ronald_B_Babu/data/SRR20707784_merge.fastq.gz'
#output_dir = '/Users/martingordon/Documents/projects/041223_Ronald_B_Babu/data'
#sample_name = 'SRR20707784'


# stage the output file path
output_path = args.output_dir + '/' + args.sample_name + '_reformat.fastq.gz'

counter=0

with gzip.open(args.fastq_path, 'rt', encoding='ascii') as fastq_file, gzip.open(output_path, 'wt') as output_file:
    for line in fastq_file:
        if line.startswith('@'):  # Header line for a read
            counter = counter + 1
            read_id = line.strip() + ' ' + args.sample_name
            seq = fastq_file.readline().strip()
            plus = fastq_file.readline().strip() 
            qual = fastq_file.readline().strip()

            output_file.write('%s\n%s\n%s\n%s\n' % (read_id, seq, plus, qual))

# write to stdout
print(f"{counter} reads detected and written to file")
