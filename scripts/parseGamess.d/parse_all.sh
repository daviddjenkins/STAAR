#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Usage: parse_all.sh storage_space_dir jobname"
    exit 3;
fi


STORAGE_DIR=$1
JOB_NAME=$2
NUM_OF_JOBS=200

i=0
while [ $i -le $NUM_OF_JOBS ]
do
    perl parseGamess.pl $STORAGE_DIR/$JOB_NAME/out/x$i $STORAGE_DIR/$JOB_NAME/STAAR/STAAR-$i.csv full > $STORAGE_DIR/$JOB_NAME/energies/STAAR-$i-wEnergy.csv 2> $STORAGE_DIR/$JOB_NAME/energies/skipped-$i.txt

    if [ -f /lustre/AQ/$JOB_NAME/inp/x$i/runs_redo.sge ]
    then
        echo "x$i has redos"; 
    fi

    if [ -n "`grep \"No results\" $STORAGE_DIR/$JOB_NAME/energies/skipped-$i.txt`" ]; then
        echo "x$i has some no results"
    fi

    i=$(( $i + 1 ))
done
