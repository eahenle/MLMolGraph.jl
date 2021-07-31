module MLMolGraph
using Distributed

using CSV, DataFrames, FIGlet, JLD2, LightGraphs, LinearAlgebra, Logging, MetaGraphs, NPZ, ProgressMeter, PyCall, Reexport, SharedArrays, SparseArrays

@reexport using PorousMaterials

import Base.show, Base.display

include("voronoi.jl")
include("vspn.jl")
include("cache_tools.jl")
include("processing.jl")
include("run_process.jl")
include("misc.jl")

function __init__()
    rc[:cache] = Dict()
    setup_cache(joinpath(pwd(), "data", "cache"))
    rc[:paths][:graphs] = joinpath(rc[:paths][:data], "graphs")
    if !isdirpath(rc[:paths][:graphs])
        mkpath(rc[:paths][:graphs])
    end
    try
        rc[:freud] = pyimport("freud")
    catch
        rc[:freud] = nothing
        @warn "Python package freud not installed"
    end
end

export  cached, xtals2primitive, bondNclassify, encode, read_targets, process_examples, clear_cache, run_process, 
        vspn_graph, voronoi_tesselation, setup_cache, VSPN_Input_Struct

end