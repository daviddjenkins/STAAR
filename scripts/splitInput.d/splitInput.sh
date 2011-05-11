#!/bin/bash
# Split the PDB list file into different lists
# so STAAR can be run in parallel on Newton
# Usage: sh splitInput.sh name list jobs

if [ $# -ne 3 ]
then
    echo "Usage: splitInput.sh name list jobs"
    exit 3;
fi

name=$1;
list=$PWD/$2;
jobs=$3;

lines=$((`wc -l $list | awk '{print $1}'`/$jobs))
echo $lines;

mkdir -p /lustre/AQ/$name/list;
cd /lustre/AQ/$name/list

split -d -l $lines $list

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

for (( i=$jobs-1; i >= 0; i--  ))
do
    mv x$i x$(($i+1))
done

ls | awk '{system("mkdir ../inp/"$1)}';
ls | awk '{system("mkdir ../out/"$1)}';
