module MLMolGraph
using Distributed

using CSV, DataFrames, FIGlet, Graphs, JLD2, LinearAlgebra, Logging, MetaGraphs, ProgressMeter, PyCall, SharedArrays, Xtals
#    , IOCapture, 
    # ,  , SparseArrays, StatsBase, 

import Base.show, Base.display

include("processing.jl")
include("run_process.jl")
include("misc.jl")
include("argument_parsing.jl")
include("export_data.jl")

function load_pydep(dep::String)
    try
        rc[Symbol(dep)] = pyimport(dep)
    catch
        rc[Symbol(dep)] = nothing
        @warn "Python dependency $dep not found."
    end
end

function __init__()
    for dir in values(rc[:paths])
        if ! isdir(dir)
            mkpath(dir)
        end
    end
    load_pydep.(["torch", "numpy", "pickle"])
end

export  run_process, parse_args, make_arg_dict

end