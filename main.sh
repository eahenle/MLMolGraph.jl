# training target
TARGET="working_capacity_temperature_swing [mmol/g]"

# types of data to produce
PROCESSING_FLAGS="--seed 42 --bonds --target_lo 0 --dataset_output _DATASET_.PKL"
# number of processes to run in data processing stage
PROCESSING_NPROC=1

# run data processing
julia -p $PROCESSING_NPROC ./main.jl $PROCESSING_FLAGS --target "$TARGET"
