#!/bin/bash
# What this does is creates the job scripts to submit to newton.
# We are doing array jobs, but since there are so many pairs,
# we can't just submit them all as a single array job.  This is
# because Newton has a limit to the largest number of array jobs,
# 10000.  So the job files are split by 10k array jobs per file.
# ALSO, Geri complained about the jobs being so short, so I needed
# combine a bunch of jobs into a single job, just running back to
# back (see the while look in the job script part).
# 
# Saying this, the script was not made to be generic. It was made
# to use once, and that's it. I split the input to STAAR into 79
# files (don't ask why I picked that number, I don't remember why.
# but I'm sure there was a reason) so we needed to make at least 
# 79 job array script files.  But because of the 10k job array 
# limitations, we ended up with something like 149 script files.
#
# Oh, and it also assumes that you split your jobs up using the
# splitInput script.  This is in the naming convention of the 
# folders.
#
# Anyway, this creates the scripts and submits the jobs

################################################################################
# CHANGE THIS IF YOU WANT TO USE THIS SCRIPT
dir=/lustre/AQ/AllPDB/
count=79
jobsperscript=10
################################################################################

# cd into the input directory
cd $dir/inp

# go through all the results files
for (( i=1; i<=$count; i++ ))
do
    ############################################################################
    # And this where it has STAAR-$i if you used a different naming scheme

    # Count the number of results
    num=$((`wc -l $dir/STAAR-$i.csv | cut -d' ' -f1` - 1))

    ############################################################################

    job=1
    start=1
    num_tasks=0
    ############################################################################
    # Again, change this if you have a different naming scheme

    # cd to the input directory for the i^th set of results
    cd x$i

    ############################################################################

    # split the runs into 10k array jobs
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

        # Output the job script
        echo -e "#$ -N Gamx$i
#$ -j y
#$ -o /dev/null
#$ -q *
#$ -pe openmpi* 1
#$ -l dedicated=4
#$ -cwd
#$ -t $start-$num_tasks:$jobsperscript

module load gamess
@ id=\$SGE_TASK_ID 
@ end=\$id
@ end += $jobsperscript

while ( \$id < \$end )
  if ( -e gamessinp-\$id.inp ) then
    mpirun -np 4 /data/AQ/bin/rungms_multi gamessinp-\$id.inp 01 1 x$i > $dir/out/x$i/gamessout-\$id.out
    rm gamessinp-\$id.F*
    rm /lustre/AQ/scr/x$i/gamessinp-\$id.*
  endif
  @ id++
end
"  > runs$i.$job.sge

        # Submit the job
        qsub runs$i.$job.sge
        job=$(($job+1))
        start=$(($num_tasks+1))

    done
    # cd back to the main input directory
    cd ../
done
