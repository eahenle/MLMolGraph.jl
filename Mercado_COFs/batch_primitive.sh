#!/bin/bash

# where to find the virtualenv for pymatgen
pymatgen_env="~/pymatgen/pymatgen_env"
# where to source the input files
inputfolder="./crystals"
# where to deposit the primitive outputs
outputfolder="./primitive_crystals"
# maximum number of concurrent processes
n_jobs=24

# get pymatgen's environment
. $pymatgen_env/bin/activate

# loop over inputs and run pymatgen to make inputs into primitive cells
joblist=($(jobs -p))
for crystal in $(ls $inputfolder)
do
    # launch jobs in background for concurrency
    python make_primitive.py $inputfolder/$crystal $outputfolder/$crystal &
    # limit the concurrent jobs to avoid system lock-up
    # not sure this works, but didn't crash, so ¯\_(ツ)_/¯
    while (( ${#joblist[*]} >= $n_jobs ))
    do
        sleep 1
        joblist=($(jobs -p))
    done
done
