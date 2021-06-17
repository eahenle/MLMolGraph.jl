"""
    `run_process(target)`
Converts inputs at `rc[:data][:crystals]` into ML model data for `target` at `rc[:data]` and `rc[:data][:graphs]`
"""
function run_process(args)
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


    # infer bonds, test quality
    good_xtals = cached("good_xtals.jld2") do
        @info "Bonding and classifying..."
        return bondNclassify(xtal_list)
    end
    @info "$(length(good_xtals)) good bonded xtals."


    # determine atom node encoding scheme
    element_to_int, max_valency = cached("encoding.jld2") do 
        @info "Determining encoding scheme..."
        return encode(good_xtals)
    end
    @info "Encoding:" element_to_int max_valency
    CSV.write("atom_to_int.csv", element_to_int)


    # get the target data
    target_df = cached("targets.jld2") do 
        @info "Reading target data..."
        return read_targets("properties.csv", good_xtals, args[:target])
    end
    CSV.write(joinpath("cofs.csv"), target_df)
    npzwrite("y.npy", target_df[:, args[:target]])


    # process the graphs into array representations
    @info "Processing examples..."
    process_examples(good_xtals, element_to_int, max_valency, args)
end
