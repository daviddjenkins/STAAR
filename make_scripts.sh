#!/bin/bash

dir=/lustre/AQ/AllPDB/
count=79

cd $dir/inp

for (( i=1; i<=$count; i++ ))
do
    num=$((`wc -l $dir/STAAR-$i.csv | cut -d' ' -f1` - 1))
    job=1
    start=1
    num_tasks=0
    cd x$i
    while [[ $num -gt 0 ]]
    do
        tmp=$num
        num=$(($num-10000))
        if [[ $num -lt 0 ]]
        then
            num_tasks=$(($num_tasks+$tmp));
        else
            num_tasks=$(($num_tasks+10000));
        fi

        echo -e "#$ -N Gamess
#$ -j y
#$ -o /dev/null
#$ -q *
#$ -pe openmpi* 1
#$ -l dedicated=4
#$ -cwd
#$ -t $start-$num_tasks

module load gamess
mpirun -np 4 /data/AQ/bin/rungms_multi gamessinp-\$SGE_TASK_ID.inp 01 1 x$i > $dir/out/x$i/gamessout-\$SGE_TASK_ID.out
rm gamessinp-\$SGE_TASK_ID.F*
rm /lustre/AQ/scr/x$i/gamessinp-\$SGE_TASK_ID.*
"  > runs$i.$job.sge

        #qsub runs$i.$job.sge
        job=$(($job+1))
        start=$(($start+num_tasks))

    done
    cd ../
done
