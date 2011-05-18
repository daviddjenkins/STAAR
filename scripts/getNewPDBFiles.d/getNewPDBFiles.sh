#!/bin/bash
# David D. Jenkins
# Script to get the new files from the PDB website
#   Usage: sh getNewPDBFiles.sh PDBList output_dir
# Very primitive, uses wget

# Check the command line args
if [[ $# -ne 2 ]]; then
    echo "Usage: sh getNewPDBFiles.sh PDBList output_dir"
    exit
fi

# Store the command line args
list=$                          # So, this is the most up to date, complete list file from www.pdb.org
outdir=$2                       # And this is where the PDB files will be stored

# This is some magical stuff.  What it does is the following:
# ls-es the PDB directory | remove the extensions | compares the list to the directory listing | takes the PDB names that are in the list, but not in the directory
#    | puts those PDB names in the URL where the PDB can be downloaded | downloads each PDB 20 at a time
ls $outdir | cut -d'.' -f1 | diff $list - | grep "<" | awk '{print "http://www.rcsb.org/pdb/files/"$2".pdb.gz"}' | xargs -P 20 -r -n 1 wget -nv -P $outdir
