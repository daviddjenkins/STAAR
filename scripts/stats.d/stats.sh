#!/bin/bash
# Gets stats on the STAAR results
# Usage: bash stats.sh staarcpp.txt STAAR.out STAARwEnergy.out

if [[ "$1" -eq "-n" ]]; then
    staarcpp='/lustre/AQ/'$2'/STAAR/staarcpp*.txt';
    STAARout='/lustre/AQ/'$2'/STAAR/STAAR*.csv';
    STAARwEnergy='/lustre/AQ/'$2'/energies/STAAR-*-wEnergy.csv';
    gamessskipped='/lustre/AQ/'$2'/energies/skipped*.txt';
else
    if [ $# -ne 4 ]
    then
        echo "Usage: bash stats.sh staarcpp.txt STAAR.out STAARwEnergy.out GAMESSskipped.txt"
        exit 3;
    fi

    # Get the command line args
    staarcpp=$1;
    STAARout=$2;
    STAARwEnergy=$3;
    gamessskipped=$4;
fi

PDBdir=/lustre/AQ/PDB

# Get total number of PDBs search again as a check
echo "= Results Summary =";
#searched=`grep ".pdb.gz" $staarcpp | grep -v "Corrected" | wc -l`;
searched=`ls --color=never $PDBdir | wc -l`
echo "Total Searched: "$searched;

echo -e "\n== Skipped PDBs ==";
# Get number of skipped due to no resolution specified
NonRes=`egrep " no resolution was specified|resolution could not be converted to a number" $staarcpp | wc -l `
echo "PDBs skipped due to no resolution specified: "$NonRes;

# Get the number of Resolution too high
ResTooHigh=`grep "resolution was too high" $staarcpp | wc -l`;
echo "PDBs skipped due to resolution too high: "$ResTooHigh;

# XRay models with multiple models
MultiModels=`grep "multiple models exist in PDB file" $staarcpp | wc -l`;
echo "Xray PDBs skipped because there were multiple models: "$MultiModels;

# Total Skipped
skipped=`egrep "Skipping.* " $staarcpp | wc -l`;
echo "Total Skipped: "$skipped;

echo -e "\n== Processed PDBs ==";
# Get the number of unique PDBs with results
PDBwResults=`grep -v "#" $STAARout | cut -d',' -f10 | sort | uniq | wc -l`;
echo "PDBs with results (pre-GAMESS): "$PDBwResults;

# PDBs with no results
PDBwoResults=$(($searched - ($skipped + $PDBwResults)));
echo "PDBs without results (pre-GAMESS): "$PDBwoResults;

echo -e "\n== Results Information - Before GAMESS ==";
# Get the number of results before GAMESS
NumResults=`grep -v "#" $STAARout | wc -l`
echo "Number of results: "$NumResults;

# Get the number of PHE-GLU results
PHEGLU=`grep -v "#" $STAARout | grep "PHE,GLU" | wc -l`;
echo "Number of PHE-GLU pairs: "$PHEGLU" | "`awk 'BEGIN{printf("%0.2f", '$PHEGLU' / '$NumResults' * 100)}'`"%";

# Get the number of PHE-ASP results
PHEASP=`grep -v "#" $STAARout | grep "PHE,ASP" | wc -l`;
echo "Number of PHE-ASP pairs: "$PHEASP" | "`awk 'BEGIN{printf("%0.2f", '$PHEASP' / '$NumResults' * 100)}'`"%";

# Get the number of PHE-PO4 results
PHEPO4=`grep -v "#" $STAARout | grep "PHE,PO4" | wc -l`;
echo "Number of PHE-PO4 pairs: "$PHEPO4" | "`awk 'BEGIN{printf("%0.2f", '$PHEPO4' / '$NumResults' * 100)}'`"%";

echo -e "\n== Results Information - After GAMESS ==";
# Get the number of results after GAMESS
NumResults=`grep -v "#" $STAARwEnergy | wc -l`
echo "Number of results: "$NumResults;

# Get the number of GAMESS skipped because of no results
NoResults=`grep "No results" $gamessskipped | wc -l`
echo "Number of pairs with no GAMESS results: "$NoResults;

# Get the number of GAMESS results that didn't converge
NoConverge=`grep "Didn't converge within 200 iterations" $gamessskipped| wc -l`
echo "Number of GAMESS runs that didn't converge: "$NoConverge;

# Get the number GAMESS results that had mix energy >.25
HighMix=`grep "abs( MIX energy )" $gamessskipped | wc -l`;
echo "Number of GAMESS results that had abs(MIX energy) > .25: "$HighMix;

# Get the number of PO4 that were rejected
PO4rej=`grep "Another PO4 result had better energy" $gamessskipped | wc -l`;
echo "Number of PO4 results that were rejected: "$(($PO4rej));

# Get the number of PHE-GLU results
PHEGLU=`grep -v "#" $STAARwEnergy | grep "PHE,GLU" | wc -l`;
echo "Number of PHE-GLU pairs: "$PHEGLU" | "`awk 'BEGIN{printf("%0.2f", '$PHEGLU' / '$NumResults' * 100)}'`"%";

# Get the number of PHE-ASP results
PHEASP=`grep -v "#" $STAARwEnergy | grep "PHE,ASP" | wc -l`;
echo "Number of PHE-ASP pairs: "$PHEASP" | "`awk 'BEGIN{printf("%0.2f", '$PHEASP' / '$NumResults' * 100)}'`"%";

# Get the number of PHE-PO4 results
PHEPO4=`grep -v "#" $STAARwEnergy | grep "PHE,PO4" | wc -l`;
echo "Number of PHE-PO4 pairs: "$PHEPO4" | "`awk 'BEGIN{printf("%0.2f", '$PHEPO4' / '$NumResults' * 100)}'`"%";
