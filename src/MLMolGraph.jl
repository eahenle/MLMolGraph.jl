module MLMolGraph
using Distributed

using CSV, DataFrames, FIGlet, JLD2, LightGraphs, LinearAlgebra, Logging, MetaGraphs, NPZ, PyCall, Reexport

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
    rc[:paths][:graphs] = joinpath(rc[:paths][:data], "Graphs")
    if !isdirpath(rc[:paths][:graphs])
        mkpath(rc[:paths][:graphs])
    end
    rc[:freud] = pyimport("freud")
end

export cached, xtals2primitive, bondNclassify, encode, read_targets, process_examples, clear_cache, run_process, vspn_graph, voronoi_tesselation

end