#!/bin/bash

source ~/.bashrc

pdbFile=${1}

mkdir -p pdbs

# read the pdbList file and print the pdb name
while read line
do
    # strip the white spaces and new lines
    pdb=$(echo $line | tr -d '\n' | tr -d ' ' | tr -d '\r' |  tr -cd '[:print:]')
    
    curl -o ./pdbs/${pdb}.pdb https://files.rcsb.org/download/${pdb}.pdb

done < ${pdbFile}
