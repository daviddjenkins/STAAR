#!/bin/bash
# Fixes failed GAMESS runs.  Sometimes the GAMESS runs
# will not recognize the input file, for some reason.
# So, we need to occasionally run this to check for failed
# runs
module load gamess

# Check the command line args
if [ $# -ne 2 ]
then
    echo "Usage: fix_failed_gamess.bash storage_space_dir jobname "
    exit 3;
fi


num_runs_per_job=10
max_time_seconds=3600 # seconds
line=1
dir=$1
scratch=$dir"/scr"
job_name=$2
inp=$dir"/"$job_name"/inp"
out=$dir"/"$job_name"/out"
tofind="This job expected the input file to be in directory"

while [[ 1 ]]
do
    queue=`qstat | awk '{print $6" "$7" "$3" "$5" "$10" "$1}' | grep -v qw | tail -n +3 | sort | head -$line | tail -1`
    d=`echo $queue | awk '{print $1" "$2}'`;
    d=`date -d"$d" +%s`
    
    curdate=`date +%s`
    elapsed=$(($curdate-$d))

    if [[ -n "$queue" ]]; then    

        if [[ $elapsed -ge $max_time_seconds ]]; then
            jobname=`echo $queue | awk '{print $3}' | sed 's/Gam//g'`
            jobfile=`echo $queue | awk '{print $5}'`

            max=$(($jobfile+$num_runs_per_job))

            cd $inp/$jobname/
            while [[ $jobfile -lt $max ]]
            do
                echo $queue
                bad=`grep "$tofind" $out/$jobname/gamessout-$jobfile.out`
            
                if [[ -n "$bad" ]]
                then
                    
                    jobfile2=$jobfile

                    while [[ $jobfile2 -lt $max ]]
                    do
                        echo $jobfile
                        mpirun -np 4 /data/AQ/bin/rungms_multi_t gamessinp-$jobfile2.inp 01 1 $jobname > $out/$jobname/gamessout-$jobfile2.out
                        rm gamessinp-$jobfile2.F*
                        rm $scratch/$jobname/gamessinp-$jobfile2.*
                        rm -r *.temp
                        jobfile2=$(($jobfile2 + 1))
                    done
                fi
                jobfile=$(($jobfile + 1))
            done
            qdel `echo $queue | awk '{print $6}'`.`echo $queue | awk '{print $5}'`
            line=$(($line + 1));
        else
            echo `date` "  :  Everything looks ok for now. Longest running: " $queue
            #break;
        fi
    else
        echo `date` "  :  Everything looks ok for now. Nothing running. "
        break;
    fi
    sleep 2m
    line=1
done
