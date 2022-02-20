#!/bin/bash
FASTA_FILE=$1
echo "Processing $FASTA_FILE file..."

# FILES=$1

# for f in $FILES/*.fa
# do
echo "(1/7) Starting CLUSTAL"
./run_clustalo.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(2/7) Starting KALIGN"
./run_kalign.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(3/7) Starting MAFFT"
./run_mafft.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(4/7) Starting MUSCLE"
./run_muscle.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(5/7) Starting MSASA BLOSUM"
./run_msasa_blosum.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(6/7) Starting MSASA GLOBAL MS"
./run_msasa_global_ms.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

echo "(7/7) Starting MSASA MATCHING"
./run_msasa_matching.sh $FASTA_FILE /home/adrian/output # &>/dev/null & disown;

# done
