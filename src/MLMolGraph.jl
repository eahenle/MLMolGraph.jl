module MLMolGraph
using Distributed

@everywhere using CSV, DataFrames, JLD2, LightGraphs, MetaGraphs, NPZ, Reexport

@reexport using Xtals # for user access to API
@everywhere using Xtals # for distributed tasks

include("cache_tools.jl")
include("processing.jl")
include("run_process.jl")

function __init__()
    rc[:cache] = Dict()
    setup_cache()
    rc[:paths][:graphs] = joinpath(rc[:paths][:data], "Graphs")
    if !isdirpath(rc[:paths][:graphs])
        mkpath(rc[:paths][:graphs])
    end
end

export cached, xtals2primitive, bondNclassify, encode, read_targets, process_examples, clear_cache, run_process

end