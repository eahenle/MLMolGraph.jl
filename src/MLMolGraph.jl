module MLMolGraph
using Distributed

using CSV, DataFrames, FIGlet, JLD2, LightGraphs, MetaGraphs, NPZ, Reexport

@reexport using PorousMaterials

include("voronoi.jl")
include("vspn.jl")
include("cache_tools.jl")
include("processing.jl")
include("run_process.jl")
include("misc.jl")

function __init__()
    rc[:cache] = Dict()
    setup_cache(joinpath(pwd(), "data", "cache"))
    rc[:paths][:graphs] = joinpath(rc[:paths][:data], "Graphs")
    if !isdirpath(rc[:paths][:graphs])
        mkpath(rc[:paths][:graphs])
    end
end

export cached, xtals2primitive, bondNclassify, encode, read_targets, process_examples, clear_cache, run_process, vspn_graph

end