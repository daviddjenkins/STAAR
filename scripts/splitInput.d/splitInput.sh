#!/bin/bash
# Split the PDB list file into different lists
# so STAAR can be run in parallel on Newton
# Note: this is assuming Newton with having the storage location /lustre/AQ
# Usage: sh splitInput.sh name list jobs project_dir

# Check the command line args
if [ $# -ne 4 ]
then
    echo "Usage: splitInput.sh name list jobs project_dir"
    exit 3;
fi

# Store the args
name=$1;
list=$PWD/$2;
jobs=$3;
dir=$4;

# Calculate the number of items that will be in each list
lines=$((`wc -l $list | awk '{print $1}'`/$jobs))
if [[ $(($lines % $jobs)) -ne 0 ]]
then
    lines=$(($lines + 1));
fi

# Make sure the list directory exists
mkdir -p $dir/$name/list;
cd $dir/$name/list

# Split up the lines
split -a 4 -d -l $lines $list

# Removes leading zeros and changes names so numbers start from 1
ls | awk -F'x' '{system( "mv "$0" x" ($2+1) )}'

# And make sure the input and output directories exist
ls | awk '{system("mkdir -p ../inp/"$1)}';
ls | awk '{system("mkdir -p ../out/"$1)}';
