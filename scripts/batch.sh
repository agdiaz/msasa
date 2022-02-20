#!/bin/bash
FASTA_FILE=$1
echo "Processing $FASTA_FILE file..."

# FILES=$1

# for f in $FILES/*.fa
# do
echo "(1/6) Starting CLUSTAL"
./run_clustalo.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(2/6) Starting KALIGN"
./run_kalign.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(3/6) Starting MAFFT"
./run_mafft.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(4/6) Starting MUSCLE"
./run_muscle.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(5/6) Starting MSASA BLOSUM"
./run_msasa_blosum.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(6/6) Starting MSASA GLOBAL MS"
./run_msasa_global_ms.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

# done
