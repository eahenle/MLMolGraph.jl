"""
    `run_process(target)`
Converts inputs at `rc[:data][:crystals]` into ML model data for `target` at `rc[:data]` and `rc[:data][:graphs]`
"""
function run_process(args)
    # get list of xtals
    @info "Finding data."
    xtal_list = cached("xtal_list.jld2") do
        return readdir(rc[:paths][:crystals])
    end
    @info "$(length(xtal_list)) inputs."

    # sample list
    if args[:samples] ≠ 0
        @info "Randomly sampling $(args[:samples])."
        xtal_list = sample(xtal_list, args[:samples], replace=false)
    end

    # process inputs to primitive cells
    cached("primitives_done.jld2") do ## TODO make a @cached macro to clean up this idiom
        if args[:primitive]
            xtals2primitive(xtal_list, args[:tolerance])
        else ## TODO make less hacky: have bondNclassify pull from different dir depending on args[:primitive]
            loadcheck(xtal_list, args[:tolerance])
        end
        return true
    end

    # infer bonds, test quality
    good_xtals = cached("good_xtals.jld2") do
        return bondNclassify(xtal_list)
    end
    @info "$(length(good_xtals)) good bonded xtals."

    # determine atom node encoding scheme
    @info "Analyzing structures."
    element_to_int, encoding_length = cached("encoding.jld2") do 
        return determine_encoding(good_xtals, args)
    end
    @info "Encoding:" element_to_int encoding_length
    CSV.write("atom_to_int.csv", element_to_int)
    npzwrite("encoding_length.npy", encoding_length)

    # get the target data
    @info "Reading target data."
    target_df = cached("targets.jld2") do 
        return read_targets("properties.csv", good_xtals, args[:target])
    end
    @info "Writing target data."
    CSV.write(joinpath("targets.csv"), target_df)
    npzwrite("y.npy", target_df[:, args[:target]])

    # process the graphs into ML inputs
    @info "Processing example data."
    process_examples(good_xtals, element_to_int, args)

    return good_xtals
end
