#!/usr/bin/env julia

# parse run flags/args
using ArgParse
argparser = ArgParseSettings()
@add_arg_table argparser begin
    "--clear_cache"
        help = "delete the cache"
        action = :store_true
    "--target", "-t"
        help = "training target"
        arg_type = String
        default = "deliverable capacity [v STP/v]"
end
args = parse_args(argparser, as_symbols=true)

using CSV, Distributed, FIGlet, NPZ
@everywhere begin
    import Pkg
    Pkg.activate(".")
end
@everywhere using MLMolGraph


FIGlet.render("MLMolGraph", FIGlet.availablefonts()[51])


if args[:clear_cache]
    @info "Clearing cache..."
    clear_cache()
end
@info "Target: $(args[:target])"


# get list of xtals
xtal_list = readdir(rc[:paths][:crystals])

@info "Found $(length(xtal_list)) inputs."


@info "Converting to primitive cells..."
xtals2primitive(xtal_list)


@info "Bonding and classifying..."
good_xtals = cached("good_xtals.jld2") do
    return bondNclassify(xtal_list)
end

@info "Kept $(length(good_xtals)) inputs."


@info "Determining encoding scheme..."
element_to_int, max_valency = cached("encoding.jld2") do 
    return encode(good_xtals)
end
CSV.write("atom_to_int.csv", element_to_int)

@info "Encoding:" element_to_int max_valency


@info "Reading target data..."
df = cached("targets.jld2") do 
    return read_targets("properties.csv", good_xtals, args[:target])
end
CSV.write(joinpath("cofs.csv"), df)
npzwrite("y.npy", df[:, args[:target]])


@info "Processing examples..."
process_examples(good_xtals, element_to_int, max_valency)

@info "Done!"
