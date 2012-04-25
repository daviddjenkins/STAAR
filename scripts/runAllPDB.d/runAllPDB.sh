#!/bin/bash
# This script does the following:
#   1. Splits the PDB list file so STAAR can run in parallel
#   2. Sets up results directory

# Check the command line args
if [ $# -ne 5 ]
then
    echo "Usage: runAllPDB.sh jobname storage_space_dir code_dir threshold resolution "
    exit 3;
fi


########################  EDIT HERE ########################

# Probably should edit these
JOB_NAME=$1
STORAGE_SPACE=$2
CODE_DIR=$3

#** STAAR parameters
THRESHOLD=$4
RESOLUTION=$5

#-----------------------------------------------------------
# Might not need to edit these
SCRIPT_HOME=`pwd`
PROJECT_DIR=$STORAGE_SPACE/$JOB_NAME
PDB_DIR=$STORAGE_SPACE/PDB
LIST_FILE=$SCRIPT_HOME/$JOB_NAME.list
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

split_up_list()
{
    rm -rf $PDB_SPLIT_LIST_DIR
    mkdir -p $PDB_SPLIT_LIST_DIR
    cd $PDB_SPLIT_LIST_DIR
    lines=$((`wc -l $LIST_FILE | awk '{print $1}'`/$NUM_STAAR_JOBS))
    if [[ $(($lines % $NUM_STAAR_JOBS)) -ne 0 ]]
    then
        lines=$(($lines + 1));
    fi

    split -a 4 -d -l $lines $LIST_FILE
    if [[ `ls $PDB_SPLIT_LIST_DIR | wc -l ` -ne $NUM_STAAR_JOBS ]]; then
        echo "Failed to split PDBs correctly!"
        echo "Num split files != Num jobs specified: " `ls $PDB_SPLIT_LIST_DIR | wc -l `" != "$NUM_STAAR_JOBS
        exit
    fi


    ls | awk -F'x' '{system( "mv "$0" x" ($2+1)) }'
    ls | awk '{system("mkdir -p ../inp/"$1)}'
    ls | awk '{system("mkdir -p ../out/"$1)}'
    cd $SCRIPT_HOME
}

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
split_up_list
echo "   ...Done"

echo "Making submission script"
create_staar_sge
echo "   ...Done"

echo "Submitting the job"
qsub $SUBMISSION_SCRIPT
echo "   ...Done"

echo "Waiting for the job to complete"
wait_for_completion
echo "   ...Done"

echo "Creating GAMESS scripts"
echo $PWD
. ../make_scripts.d/make_scripts.sh $JOB_NAME $STORAGE_SPACE $NUM_STAAR_JOBS
cd $SCRIPT_HOME
echo "   ...Done"
