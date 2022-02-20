#!/bin/bash

# ./run_muscle.sh ~/data/SELECTED/small_short.fa ~/output

EXECUTIONS_PER_TOOL=1
INPUT_FASTA=$1
OUTPUT_FOLDER=$2

filename=$(basename -- "$INPUT_FASTA")
extension="${filename##*.}"
filename="${filename%.*}"

mkdir -p $OUTPUT_FOLDER/$filename
LOG_FILENAME=$OUTPUT_FOLDER/$filename/$filename.muscle.log

for ((i=1;i<=EXECUTIONS_PER_TOOL;i++));
do
	echo "MUSCLE Execution # $i started"
	# mkdir -p $OUTPUT_FOLDER/$filename
	muscle -quiet -in $INPUT_FASTA -out $OUTPUT_FOLDER/$filename/$i.muscle.fa -log $LOG_FILENAME
done

cat $LOG_FILENAME

echo "DONE"
