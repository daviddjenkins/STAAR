#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: get_redos.sh /path/to/results/folder"
    exit 3;
fi

script_dir=`pwd`
cd $1;

for file in `ls *.redo`
do
    perl $script_dir/get_energies.pl $file
done

cd - > /dev/null
