# set up environment and paths
using CSV, Distributed, NPZ
@everywhere begin
    import Pkg
    Pkg.activate(".")
end
@everywhere using COFProcessing


@info "\n\n\tCOF Processing\n\n\n"


if length(ARGS) == 1 && ARGS[1] == "reset"
    @info "Clearing cache..."
    clear_cache()
end


# choose target
target_symbol = Symbol("deliverable capacity [v STP/v]")

@info "Target: $target_symbol"


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
    return read_targets("properties.csv", good_xtals, target_symbol)
end
CSV.write(joinpath("cofs.csv"), df)
npzwrite("y.npy", df[:, target_symbol])


@info "Processing examples..."
process_examples(good_xtals, element_to_int, max_valency)

@info "Done!"
