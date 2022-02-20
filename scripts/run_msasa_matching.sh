#!/bin/bash

# ./run_msasa.sh ~/data/SELECTED/small_short.fa ~/output

EXECUTIONS_PER_TOOL=1
INPUT_FASTA=$1
OUTPUT_FOLDER=$2
MSASA=/home/adrian/msasa/src/msa.py

filename=$(basename -- "$INPUT_FASTA")
extension="${filename##*.}"
filename="${filename%.*}"

mkdir -p $OUTPUT_FOLDER/$filename
LOG_FILENAME=$OUTPUT_FOLDER/$filename/$filename.matching.log

echo "file;initial_energy;final_energy" > $LOG_FILENAME

for ((i=1;i<=EXECUTIONS_PER_TOOL;i++));
do
	OUTPUT_MSA_FILENAME=$OUTPUT_FOLDER/$filename/$i.msasa_matching.fa
	PLOT_BEST_FILENAME=$OUTPUT_FOLDER/$filename/$i.msasa_matching_best.png
	PLOT_TEMP_FILENAME=$OUTPUT_FOLDER/$filename/$i.msasa_matching_temp.png
	echo "MSASA Execution # $i started"
	python3 $MSASA --input $INPUT_FASTA --output $OUTPUT_MSA_FILENAME --comparer matching --n-iterations 10000 --output-best-plot $PLOT_BEST_FILENAME --output-temp-plot $PLOT_TEMP_FILENAME --optimization min --temperature 1 >> $LOG_FILENAME
done

LOG_FILENAME_SORTED=$LOG_FILENAME.sorted

cat $LOG_FILENAME | tail -n +2 | sort -r --field-separator=";" --key=3 > $LOG_FILENAME_SORTED
cat $LOG_FILENAME_SORTED


BEST=$(cat $LOG_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )

ln --symbolic --force $BEST $OUTPUT_FOLDER/$filename/best.msasa_matching.fa
echo "DONE"
