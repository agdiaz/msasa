#!/bin/bash

# ./run_kalign.sh ~/data/SELECTED/small_short.fa ~/output

EXECUTIONS_PER_TOOL=1
INPUT_FASTA=$1
OUTPUT_FOLDER=$2
MSASA=/home/adrian/msasa/src/msa.py

filename=$(basename -- "$INPUT_FASTA")
extension="${filename##*.}"
filename="${filename%.*}"

mkdir -p $OUTPUT_FOLDER/$filename
LOG_FILENAME=$OUTPUT_FOLDER/$filename/$filename.kalign.log

for ((i=1;i<=EXECUTIONS_PER_TOOL;i++));
do
	echo "KAlign Execution # $i started"
	kalign -quiet -input $INPUT_FASTA -output $OUTPUT_FOLDER/$filename/$i.kalign.fa -format fasta > LOG_FILENAME
done

echo "DONE"
