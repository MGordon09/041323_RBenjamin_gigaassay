#!/bin/bash

fastq=$1

export header=$(gzip -cd $1 | head -n 1)
export sm=$(echo $header | head -n 1 | cut -f 2 -d "."| cut -f 1 -d " ") #| sed 's/@//' | sed 's/:/_/g')
export id=$(echo $header | head -n 1 | cut -f 1 -d " "| cut -f 1 -d "." | sed 's/@//' )
#echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tPL:ILLUMINA"