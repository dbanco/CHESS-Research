#!/bin/bash

dataPath="/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/"
outPath="nfs/chess/user/dbanco/outputs/"
params="params_string"  # Make sure this is formatted correctly to be parsed in the Python script

# Initial values
t=scan1
spotInds=(0 1 2 3)  # List of spot indices
fname1="path_to_fname1"
fname2="path_to_fname2"

# Submit jobs for each spot
for k in ${!spotInds[@]}; do
    s=${spotInds[$k]}
    qsub -N processSpotJob_$k -o output_$k.txt -e error_$k.txt -b y "python process_spot.py $k $s $t $params $outPath $fname1 $fname2"
done
