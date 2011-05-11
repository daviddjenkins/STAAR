#!/bin/bash

for line in $(< $1);do
    if [ ! -f "/lustre/AQ/PDB/$line.pdb.gz" ]; then
        echo "http://www.rcsb.org/pdb/files/$line.pdb.gz";
    fi
done
