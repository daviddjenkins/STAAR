#!/bin/bash
# David D. Jenkins
# Script to get the new files from the PDB website
#   Usage: sh getNewPDBFiles.sh PDBList output_dir
# Very primitive, uses wget

# Check the command line args
if [[ $# -ne 1 ]]; then
    echo "Usage: sh getNewPDBFiles.sh output_dir"
    exit
fi

# Store the command line args
outdir=$1                       # And this is where the PDB files will be stored
list=PDBList`date "+%Y%m%d"`.txt;

# Download the latest listing
wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt -O - -q | cut -c 1-4 | sort | tr '[:lower:]' '[:upper:]' > $list.whole

# Get the differences between what we already have and the new list
ls --color=never $outdir | cut -d'.' -f1 | sort | diff $list.whole - | grep "<" > $list.new

# And get the depreciated PDBs
ls --color=never $outdir | cut -d'.' -f1 | sort | diff $list.whole - | grep ">" > $list.del

# This is some magical stuff.  What it does is the following:
# ls-es the PDB directory | remove the extensions | compares the list to the directory listing | takes the PDB names that are in the list, but not in the directory
#    | puts those PDB names in the URL where the PDB can be downloaded | downloads each PDB 20 at a time
#ls $outdir | cut -d'.' -f1 | diff $list - | grep "<" | awk '{print "http://www.rcsb.org/pdb/files/"$2".pdb.gz"}' | xargs -P 20 -r -n 1 wget -nv -P $outdir

# But now, we don't need all of those.  we just stored the list to the disk for use later
awk '{print "http://www.rcsb.org/pdb/files/"$2".pdb.gz"}' $list.new | xargs -P 20 -r -n 1 wget -nv -P $outdir.update

# Delete the extra depreciated files
awk '{print $outdir/$2}' $list.del | xargs -P 20 -r -n 1 rm 

