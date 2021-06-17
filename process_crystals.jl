#!/usr/bin/env julia
using ArgParse
argparser = ArgParseSettings(
    prog="process_crystals", 
    description="reads .cif inputs, converts to primitive cell, infers bonds, and converts valid bonded structures to ML inputs for GNNs"
    )
@add_arg_table argparser begin
    "--target", "-t"
        help = "training target"
        arg_type = Symbol
        default = Symbol("deliverable capacity [v STP/v]")
    "--data", "-d"
        help = "root of data tree. inputs loaded from `crystals` subdirectory."
        arg_type = String
        default = joinpath(pwd(), "data")
    "--bonds", "-b"
        help = "write ML inputs for the bonding graph"
        action = :store_true
    "--angles", "-a"
        help = "write ML inputs for bonding angles"
        action = :store_true
    "--vspn", "-v"
        help = "write ML inputs for Voronoi sphere pore network"
        action = :store_true
    "--forcefield", "-f"
        help = "the forcefield to use for VSPN"
        arg_type = String
        default = ""
    "--probe", "-p"
        help = "the molecule to use as the probe in VSPN"
        arg_type = String
        default = ""
    "--clear_cache", "-c"
        help = "delete the cache"
        action = :store_true
end
args = parse_args(argparser, as_symbols=true)

using CSV, FIGlet, MLMolGraph, NPZ

FIGlet.render("MLMolGraph", FIGlet.availablefonts()[64])
@info "Loading data from $(args[:data])"
@info "Cache at $(rc[:paths][:data])/cache"
@info "Target: $(args[:target])"

set_paths(args[:data])
if args[:clear_cache]
    @info "Clearing cache..."
    clear_cache()
end

run_process(args)

@info "Done!"
