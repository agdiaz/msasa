#!/bin/bash

FILES=$1

for f in $FILES/*.fa
do
	echo "Processing $f file..."
	./run_msasa.sh $f /home/adrian/output &>/dev/null & disown;
	./run_clustalo.sh $f /home/adrian/output &>/dev/null & disown;
	./run_kaling.sh $f /home/adrian/output &>/dev/null & disown;
	./run_mafft.sh $f /home/adrian/output &>/dev/null & disown;
done
