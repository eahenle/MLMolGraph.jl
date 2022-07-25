# training target
TARGET="working_capacity_vacuum_swing [mmol/g]"

# types of data to produce
PROCESSING_FLAGS="--bonds"
# number of processes to run in data processing stage
PROCESSING_NPROC=2

# run data processing
julia -p $PROCESSING_NPROC ./main.jl $PROCESSING_FLAGS --target "$TARGET"
