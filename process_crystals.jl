#!/usr/bin/env julia
using ArgParse
argparser = ArgParseSettings(
    prog="process_crystals", 
    description="reads .cif inputs, converts to primitive cell, infers bonds, and converts valid bonded structures to ML inputs for GNNs"
    )
@add_arg_table argparser begin
    "--clear_cache", "-c"
        help = "delete the cache"
        action = :store_true
    "--target", "-t"
        help = "training target"
        arg_type = String
        default = "deliverable capacity [v STP/v]"
    "--data", "-d"
        help = "root of data tree. inputs loaded from `crystals` subdirectory."
        arg_type = String
        default = joinpath(pwd(), "data")
end
args = parse_args(argparser, as_symbols=true)

using CSV, FIGlet, MLMolGraph, NPZ

FIGlet.render("MLMolGraph", FIGlet.availablefonts()[51])
@info "Loading data from $(args[:data])"
@info "Cache at $(rc[:paths][:data]/cache)"
@info "Target: $(args[:target])"

set_paths(args[:data])
if args[:clear_cache]
    @info "Clearing cache..."
    clear_cache()
end

run_process(args[:target])

@info "Done!"
