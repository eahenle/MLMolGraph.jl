#!/bin/bash

# get number of procs from CLI; if not specified, default to single-threaded
if [[ $1 == "" ]]
then
    NB_PROC_ARG=""
else
    NB_PROC_ARG="-p $1"
fi

echo "julia $NB_PROC_ARG process_crystals.jl args..."

# run the processing script w/ provided input flags
julia $NB_PROC_ARG process_crystals.jl --samples 10000 --bonds --vspn --forcefield UFF --probe He --clear_cache
