#!/bin/bash
# Split the PDB list file into different lists
# so STAAR can be run in parallel on Newton
# Note: this is assuming Newton with having the storage location /lustre/AQ
# Usage: sh splitInput.sh name list jobs

# Check the command line args
if [ $# -ne 3 ]
then
    echo "Usage: splitInput.sh name list jobs"
    exit 3;
fi
##################################################
# Change if you want to store elsewhere
dir="/lustre/AQ/"
##################################################

# Store the args
name=$1;
list=$PWD/$2;
jobs=$3;

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
split -d -l $lines $list

# And here I got lazy.  This just removes the leading zeros
# will not work if you have more than 9999 splits
if [ $jobs -ge 10 ]
then
    r=`ls x0* | sed 's/x0//g' | awk '{system("mv x0"$1" x"$1)}'`;
elif [ $jobs -ge 100 ]
then
    r=`ls x00* | sed 's/x00//g' | awk '{system("mv x00"$1" x"$1)}'`;
    r=`ls x0* | sed 's/x0//g' | awk '{system("mv x0"$1" x"$1)}'`;
elif [ $jobs -ge 1000 ]
then
    r=`ls x000* | sed 's/x000//g' | awk '{system("mv x000"$1" x"$1)}'`;
    r=`ls x00* | sed 's/x00//g' | awk '{system("mv x00"$1" x"$1)}'`;
    r=`ls x0* | sed 's/x0//g' | awk '{system("mv x0"$1" x"$1)}'`;
fi

# Again, got lazy...
# make the lists be 1 indexed instead of 0
for (( i=$jobs-1; i >= 0; i--  ))
do
    mv x$i x$(($i+1))
done

# And make sure the input and output directories exist
ls | awk '{system("mkdir -p ../inp/"$1)}';
ls | awk '{system("mkdir -p ../out/"$1)}';
