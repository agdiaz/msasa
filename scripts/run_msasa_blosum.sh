#!/bin/bash

# ./run_msasa.sh ~/data/SELECTED/small_short.fa ~/output/small_short.fa.out1

EXECUTIONS_PER_TOOL=30
INPUT_FASTA=$1
OUTPUT_FOLDER=$2
MSASA=/home/adrian/msasa/src/msa.py

filename=$(basename -- "$INPUT_FASTA")
extension="${filename##*.}"
filename="${filename%.*}"

mkdir -p $OUTPUT_FOLDER/$filename
LOG_FILENAME=$OUTPUT_FOLDER/$filename/$filename.log
# echo "Log file with global results: $LOG_FILENAME"

echo "file;initial_energy;final_energy" > $LOG_FILENAME
# date >> $LOG_FILENAME

for i in {1..5}
do
	OUTPUT_MSA_FILENAME=$OUTPUT_FOLDER/$filename/$i.msasa-blosum.fa
	PLOT_FILENAME=$OUTPUT_FOLDER/$filename/$i.blosum.png
	echo "MSASA Execution # $i started"
	# mkdir -p $OUTPUT_FOLDER/$filename
	python3 $MSASA $INPUT_FASTA $OUTPUT_MSA_FILENAME blosum 5000 $PLOT_FILENAME >> $LOG_FILENAME
done

LOG_FILENAME_SORTED=$LOG_FILENAME.sorted

cat $LOG_FILENAME | tail -n +2 | sort -r --field-separator=";" --key=3 > $LOG_FILENAME_SORTED
cat $LOG_FILENAME_SORTED


BEST=$(cat $LOG_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )

ln -s $BEST $OUTPUT_FOLDER/$filename/best.msasa.fa
echo "DONE"
