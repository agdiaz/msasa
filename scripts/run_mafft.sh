#!/bin/bash

# ./run_mafft.sh ~/data/SELECTED/small_short.fa ~/output

EXECUTIONS_PER_TOOL=1
INPUT_FASTA=$1
OUTPUT_FOLDER=$2

filename=$(basename -- "$INPUT_FASTA")
extension="${filename##*.}"
filename="${filename%.*}"

mkdir -p $OUTPUT_FOLDER/$filename

for ((i=1;i<=EXECUTIONS_PER_TOOL;i++));
do
	echo "MAFFT Execution # $i started"
	mafft --auto --quiet $INPUT_FASTA > $OUTPUT_FOLDER/$filename/$i.mafft.fa
done

echo "DONE"
