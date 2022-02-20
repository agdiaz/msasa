#!/bin/bash

# ./run_clustalo.sh ~/data/SELECTED/small_short.fa ~/output

EXECUTIONS_PER_TOOL=1
INPUT_FASTA=$1
OUTPUT_FOLDER=$2

filename=$(basename -- "$INPUT_FASTA")
extension="${filename##*.}"
filename="${filename%.*}"

mkdir -p $OUTPUT_FOLDER/$filename
LOG_FILENAME=$OUTPUT_FOLDER/$filename/$filename.clustalo.log

for ((i=1;i<=EXECUTIONS_PER_TOOL;i++));
do
	echo "ClustalOmega Execution # $i started"
	clustalo --force -i $INPUT_FASTA -o $OUTPUT_FOLDER/$filename/$i.clustalo.fa > $LOG_FILENAME
done

cat $LOG_FILENAME

echo "DONE"
