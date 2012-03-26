#!/bin/bash
# This script does the following:
#   1. Splits the PDB list file so STAAR can run in parallel
#   2. Sets up results directory

########################  EDIT HERE ########################

# Probably should edit these
JOB_NAME=PDB20120325
STORAGE_SPACE=/lustre/AQ
CODE_DIR=/home/djenki11/scratch/repo/

#** STAAR parameters
THRESHOLD=7
RESOLUTION=3

#-----------------------------------------------------------
# Might not need to edit these
PROJECT_DIR=$STORAGE_SPACE/$JOB_NAME
PDB_DIR=$STORAGE_SPACE/PDB
LIST_FILE=$JOB_NAME.list
NUM_STAAR_JOBS=200
PDB_SPLIT_LIST_DIR=$PROJECT_DIR/list
STAAR_OUT_DIR=$PROJECT_DIR/STAAR
SUBMISSION_SCRIPT=$PROJECT_DIR/STAAR.sge
GAMESS_INP_DIR=$PROJECT_DIR/inp
GAMESS_OUT_DIR=$PROJECT_DIR/out

# May not need this if you want to supply your own list
ls --color=never $PDB_DIR | awk -F "." '{print $1}' > $LIST_FILE
#-----------------------------------------------------------

############################################################

create_staar_sge()
{
    echo "#!/bin/bash

# Name of the job
#\$ -N $JOB_NAME

# Pass all of the enivornment variables to job
#\$ -V

# Set the file to which STDERR will be written
#\$ -e $STAAR_OUT_DIR/STAAR.err

# Set the file to which STDOUT will be written
#\$ -o $STAAR_OUT_DIR/STAAR.out

# Set what queue we are going to run it in
#\$ -q medium*

# Move to the current working directory
#\$ -cwd

# Email address to send completion notice to
#\$ -M david.d.jenkins@gmail.com

# Send an email on suspension of the job
#\$ -m bes

# Set up an array job
#\$ -t 1-$NUM_STAAR_JOBS

cd $CODE_DIR

$CODE_DIR/staar --residues \"PHE;GLU,ASP\" --threshold $THRESHOLD --resolution $RESOLUTION --out $STAAR_OUT_DIR/STAAR-\$SGE_TASK_ID.csv --pdbdir $PDB_DIR --pdblist $PDB_SPLIT_LIST_DIR/x\$SGE_TASK_ID --gamess $GAMESS_INP_DIR/x\$SGE_TASK_ID > $STAAR_OUT_DIR/staarcpp-\$SGE_TASK_ID.txt

" > $SUBMISSION_SCRIPT
}

wait_for_completion()
{
    until [ -z "`qstat | grep "PDB2"`" ]; do
        sleep 3m
    done
}

echo "Making directory structure..."
mkdir -p $STAAR_OUT_DIR
mkdir -p $PDB_SPLIT_LIST_DIR
mkdir -p $GAMESS_INP_DIR
mkdir -p $GAMESS_OUT_DIR
echo "   ...Done"

echo "Splitting PDB files into different lists and setting up job dir..."
. ../splitInput.d/splitInput.sh $JOB_NAME $LIST_FILE $NUM_STAAR_JOBS $STORAGE_SPACE
echo "   ...Done"

echo "Making submission script"
create_staar_sge
echo "   ...Done"

echo "Submitting the job"
#qsub $SUBMISSION_SCRIPT
echo "   ...Done"

echo "Waiting for the job to complete"
wait_for_completion
echo "   ...Done"

