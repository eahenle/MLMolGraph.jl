"""
    `run_process(target)`
Converts inputs at `rc[:data][:crystals]` into ML model data for `target` at `rc[:data]` and `rc[:data][:graphs]`
"""
function run_process(args)
    # check args
    valid, msg = validate_args(args)
    if !valid
        error(msg)
    end

    # get list of xtals
    xtal_list = cached("xtal_list.jld2") do
        return readdir(rc[:paths][:crystals])
    end
    @info "$(length(xtal_list)) inputs."


    # process inputs to primitive cells
    cached("primitives_done.jld2") do
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
    element_to_int, max_valency, encoding_length = cached("encoding.jld2") do 
        return encode(good_xtals)
    end
    @info "Encoding:" element_to_int max_valency encoding_length
    CSV.write("atom_to_int.csv", element_to_int)
    npzwrite("encoding_length.npy", encoding_length)


    # get the target data
    target_df = cached("targets.jld2") do 
        @info "Reading target data..."
        return read_targets("properties.csv", good_xtals, args[:target])
    end
    CSV.write(joinpath("targets.csv"), target_df)
    npzwrite("y.npy", target_df[:, args[:target]])


    # process the graphs into ML inputs
    process_examples(good_xtals, element_to_int, max_valency, args)

    if args[:env] == "julia" && args[:vspn]
        @info "Preparing VSPN data dictionary..."
        # pack data into array for VSPN dataloader
        cached("vspns.jld2") do
            vspns = VSPN_Input_Struct[]
            for (i, xtal) ∈ enumerate(good_xtals)
                xtal_name = split(xtal, ".cif")[1]
                @load joinpath(rc[:cache][:vspn], "$xtal_name.jld2") obj
                g, X = obj
                y = target_df[i, args[:target]]
                push!(vspns, VSPN_Input_Struct(
                    xtal_name,
                    adjacency_matrix(g),
                    sparse(X),
                    y
                ))
            end
            return vspns
        end
    end
end
