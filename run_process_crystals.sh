#!/bin/bash

# get number of threads from CLI; if not specified, default to 2
if [[ $1 == "" ]]
then
    NB_PROC_ARG=""
else
    NB_PROC_ARG="-p $1"
fi

echo "julia $NB_PROC_ARG process_crystals.jl args..."

# run the processing script w/ provided input flags
julia $NB_PROC_ARG process_crystals.jl --bonds --vspn --forcefield UFF --probe CH4 --sample 2000 --clear_cache --verbose

## ? what's the deal w/ warning about requiring OffsetArrays? Isn't resolved by `]add`
