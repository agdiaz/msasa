#!/bin/bash

DIRECTORY=$1

for fasta_file in $DIRECTORY/*tfa; do 
	echo "Processing: $fasta_file";
	cat $fasta_file | grep '>' | wc -l > $fasta_file.stats.txt.count.txt
done
