#!/usr/bin/env julia


using CSV, FIGlet, MLMolGraph, NPZ


# command-line options/help
using ArgParse
argparser = ArgParseSettings(
    prog="process_crystals", 
    description="reads .cif inputs, converts to primitive cell, infers bonds, and converts valid bonded structures to ML inputs for GNNs"
    )
@add_arg_table argparser begin
    "--clear-cache"
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


# banner
FIGlet.render("MLMolGraph", FIGlet.availablefonts()[51])


# path and cache setup
set_paths(args[:data])
if args[:clear_cache]
    @info "Clearing cache..."
    clear_cache()
end


# parameter report
@info "Loading data from $(args[:data])"
@info "Cache at $(rc[:paths][:data]/cache)"
@info "Target: $(args[:target])"


# get list of xtals
xtal_list = cached("xtal_list.jld2") do
    readdir(rc[:paths][:crystals])
end
@info "$(length(xtal_list)) inputs."


# process inputs to primitive cells
cached("primitives_done.jld2") do
    @info "Converting to primitive cells..."
    xtals2primitive(xtal_list)
    return true
end


good_xtals = cached("good_xtals.jld2") do
    @info "Bonding and classifying..."
    return bondNclassify(xtal_list)
end
@info "$(length(good_xtals)) good bonded xtals."


element_to_int, max_valency = cached("encoding.jld2") do 
    @info "Determining encoding scheme..."
    return encode(good_xtals)
end
@info "Encoding:" element_to_int max_valency
CSV.write("atom_to_int.csv", element_to_int)


target_df = cached("targets.jld2") do 
    @info "Reading target data..."
    return read_targets("properties.csv", good_xtals, args[:target])
end
CSV.write(joinpath("cofs.csv"), target_df)
npzwrite("y.npy", target_df[:, args[:target]])


@info "Processing examples..."
process_examples(good_xtals, element_to_int, max_valency)


@info "Done!"
