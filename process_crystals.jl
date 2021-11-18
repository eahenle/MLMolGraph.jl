#!/usr/bin/env julia

using ArgParse, Distributed

argparser = ArgParseSettings(
    prog="process_crystals", 
    description="reads .cif inputs, converts to primitive cell, infers bonds, and converts valid bonded structures to ML inputs for GNNs"
)

@add_arg_table argparser begin
    "--angles", "-a"
        help = "write ML inputs for bonding angles"
        action = :store_true

    "--bonds", "-b"
        help = "write ML inputs for the bonding graph"
        action = :store_true

    "--cache", "-c"
        help = "cache location"
        arg_type = String
        default = joinpath(pwd(), "data", "cache")

    "--data", "-d"
        help = "root of data tree. inputs loaded from `crystals` subdirectory."
        arg_type = String
        default = joinpath(pwd(), "data")
    
    "--verbose", "-e"
        help = "whether or not to print extra @info output"
        action = :store_true

    "--forcefield", "-f"
        help = "the forcefield to use for VSPN"
        arg_type = String
        default = ""

    "--vectors", "-g"
        help = "write bond vectors"
        action = :store_true

    "--tolerance", "-n"
        help = "set net charge tolerance for Crystal()"
        arg_type = Float64
        default = 1e-4

    "--one_hot", "-o"
        help = "choose one-hot encodings"
        arg_type = String
        default = "n" # for "normal"
         
    "--probe", "-p"
        help = "the molecule to use as the probe in VSPN"
        arg_type = String
        default = ""
       
    "--primitive", "-r"
        help = "reduce input to primitive cell"
        arg_type = Bool
        default = true

    "--samples", "-s"
        help = "number of inputs to sample (0 -> use all)"
        arg_type = Int
        default = 0

    "--target", "-t"
        help = "training target"
        arg_type = Symbol
        default = Symbol("working_capacity_vacuum_swing [mmol/g]")

    "--vspn", "-v"
        help = "write ML inputs for Voronoi sphere pore network"
        action = :store_true

    "--clear_cache", "-x"
        help = "delete the cache"
        action = :store_true
end

args = parse_args(argparser, as_symbols=true)

@everywhere begin ## pre-registration dev hack. makes wall of annoying console output.
    import Pkg
    Pkg.activate(".")
end

@everywhere using MLMolGraph

MLMolGraph.banner()

# check args
valid, msg = MLMolGraph.validate_args(args)
if !valid
    error(msg)
end

@info "Program parameters" args...

set_paths(args[:data])
setup_cache(args[:cache])
if args[:clear_cache]
    @info "Clearing cache."
    clear_cache()
end

run_process(args)

@info "Done!"
