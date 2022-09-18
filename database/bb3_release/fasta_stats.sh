#!/bin/bash

DIRECTORY=$1

for fasta_file in $DIRECTORY/*.tfa; do 
	echo "Processing $fasta_file";
	cat $fasta_file | awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' > $fasta_file.stats
done
