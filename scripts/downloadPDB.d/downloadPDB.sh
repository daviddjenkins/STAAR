#! /bin/bash
# Downloads the PDBs from a list if they do not already exist

if [[ $# -ne 2 ]]; then
    echo "Usage: sh getNewPDBFiles.sh PDBList output_dir"
    exit
fi

# This is an option to get compressed or uncompressed files
# hidden option because it isn't necessary for most purposes
comp=1

# This is the address that the PDB files will be retrieved from
site="http://www.rcsb.org/pdb/files/"

# You shouldn't have to edit anything below this line (hopefully)

# List of PDB files
pdblist=$1

# output directory
outdir=$2

# Set the extension based on compression option
if [[ "$comp" -eq "1" ]]; then
    ext="pdb.gz"
else
    ext="pdb"
fi

echo "Latest PDB files: "
# Now, we go through each line of the PDB list file
for line in $(< $pdblist);do
    # If it doesn't exists already, we need to download it
    if [ ! -f "$outdir/$line.$ext" ]; then
        wget -P $outdir $site/$line.$ext
        echo $line.$ext
    fi
done

