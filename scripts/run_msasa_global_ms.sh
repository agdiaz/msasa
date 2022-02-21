#!/bin/bash

# ./run_msasa.sh ~/data/SELECTED/small_short.fa ~/output

EXECUTIONS_PER_TOOL=30
INPUT_FASTA=$1
OUTPUT_FOLDER=$2
MSASA=/home/adrian/msasa/src/msa.py

filename=$(basename -- "$INPUT_FASTA")
extension="${filename##*.}"
filename="${filename%.*}"

mkdir -p $OUTPUT_FOLDER/$filename
LOG_FILENAME=$OUTPUT_FOLDER/$filename/$filename.global_ms.log

echo "file;initial_energy;final_energy" > $LOG_FILENAME

for ((i=1;i<=EXECUTIONS_PER_TOOL;i++));
do
	OUTPUT_MSA_FILENAME=$OUTPUT_FOLDER/$filename/$i.msasa_global_ms.fa
	PLOT_BEST_FILENAME=$OUTPUT_FOLDER/$filename/$i.msasa_global_ms_best.png
	PLOT_TEMP_FILENAME=$OUTPUT_FOLDER/$filename/$i.msasa_global_ms_temp.png
	echo "MSASA Execution # $i started"
	python3 $MSASA --input $INPUT_FASTA --output $OUTPUT_MSA_FILENAME --comparer global_ms_min --n-iterations 50 --output-best-plot $PLOT_BEST_FILENAME --output-temp-plot $PLOT_TEMP_FILENAME --optimization min >> $LOG_FILENAME
done

LOG_FILENAME_SORTED=$LOG_FILENAME.sorted

cat $LOG_FILENAME | tail -n +2 | sort -r --field-separator=";" --key=3 > $LOG_FILENAME_SORTED
cat $LOG_FILENAME_SORTED


BEST=$(cat $LOG_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )

ln --symbolic --force $BEST $OUTPUT_FOLDER/$filename/best.msasa_global_ms.fa
echo "DONE"
