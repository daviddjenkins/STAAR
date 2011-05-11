#! /bin/bash
# Parses the PDB Select list from 
#   http://bioinfo.tg.fh-giessen.de/pdbselect/recent.pdb_select25
# if no file is given
# Puts them in the form of
#   PDB_ID Chain_ID_1,Chain_ID_2
#
# Usage: 
#    parsePDBSelect.sh 
#  or
#    parsePDBSelect.sh /path/to/PDBSelect/file

URL="http://bioinfo.tg.fh-giessen.de/pdbselect/recent.pdb_select25";

if [ $# -eq 0 ]
then
    curl -O $URL;
    file="recent.pdb_select25";
else
    file=$1;
fi

# Get the PDB ids
grep -v '#' $file | cut -c 9-12 > temp1.txt; 

# Get the chains
grep -v '#' $file | cut -c 13 > temp2.txt; 

# Combine and sort them.  Then put chains that are from
# the same PDB in the same line separated by a comma
cat temp2.txt | paste temp1.txt - | awk '
   BEGIN {FS=OFS="\t"}
   !A[$1] {A[$1] = $0; next}
   {A[$1] = A[$1] "," $2}
   END {for(i in A) {print A[i]}
}';

# Remove extra files
rm temp{1,2}.txt;

if [ $# -eq 0 ]
then
    rm $file;
fi

