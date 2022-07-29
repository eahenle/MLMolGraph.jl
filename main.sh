# training target
TARGET=""

# types of data to produce
PROCESSING_FLAGS="--seed 42 --bonds --dataset_output _DATASET_.PKL"
# number of processes to run in data processing stage
PROCESSING_NPROC=1

# run data processing
julia -p $PROCESSING_NPROC ./main.jl $PROCESSING_FLAGS --target "$TARGET"
