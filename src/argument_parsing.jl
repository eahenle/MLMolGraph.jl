using ArgParse, Random

# create the argument parser object
argparser = ArgParseSettings(
    prog="MLMolGraph Data Processor", # program name to display
    description="reads .cif inputs, converts to primitive cell, infers bonds, and converts valid bonded structures to ML inputs"
)

# build the argument table
include("argument_table.jl")


# simple function for parsing command line args
import ArgParse.parse_args
parse_args() = parse_args(argparser, as_symbols=true)


"""
    args::Dict{Symbol, Any} = make_arg_dict(target::String, args::Dict{Symbol, Any})
    args::Dict{Symbol, Any} = make_arg_dict(args::Dict{Symbol, Any})
"""
function make_arg_dict(target::String=DEFAULT_ARGS[:target], args::Dict{Symbol, Any}=Dict{Symbol, Any}())::Dict{Symbol, Any}
    args = merge(DEFAULT_ARGS, args)
    args[:target] = target
    return args
end

make_arg_dict(args::Dict{Symbol, Any}) = make_arg_dict(args[:target], args)
