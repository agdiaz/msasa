#!/bin/bash

BASE_DIR=$PWD
DIRECTORY=$BASE_DIR/$1
for f in $DIRECTORY/*.tfa; do
    file=${f##*/}
    base=${file%.*}
    echo "Base file: $DIRECTORY/$base"
    export INPUT_SEQUENCES=$DIRECTORY/$base
    export ITERATIONS=30
    make run-workflow
done
