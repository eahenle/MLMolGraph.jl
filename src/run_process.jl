"""
    `run_process(args)`
Converts inputs at `rc[:data][:crystals]` into ML model data at `rc[:data]`
"""
function run_process(args; temp_dir=nothing)
    # create temporary working directory
    temp_dir = isnothing(temp_dir) ? tempdir() : temp_dir
    mkpath(temp_dir)
    @assert isdir(temp_dir) "Failed to create temporary directory at $temp_dir"

    # copy selected rows into smaller CSV file (for convenience)
    @info "Reading target data."
    target_csv = joinpath(rc[:paths][:data], "properties.csv")
    df = CSV.read(target_csv, DataFrame)
    @assert nrow(df) > 0 "$target_csv contains no examples."
    @assert "name" in names(df) "$target_csv has no 'name' column."

    # get list of xtals
    xtal_dir = rc[:paths][:crystals]
    all_files = readdir(xtal_dir)
    file_is_cif = @showprogress "Finding example data " pmap(x -> contains(x, ".cif"), all_files)
    validation_progress = Progress(6, "Validating example data ")
    # filter for just the CIF data files
    cif_files = all_files[file_is_cif .== 1]
    next!(validation_progress)
    @assert length(cif_files) != 0 "$(xtal_dir) contains no CIF files."
    next!(validation_progress)
    mof_names = [String(split(cif_file, ".cif")[1]) for cif_file in cif_files]
    next!(validation_progress)
    
    found_cif = intersect(mof_names, df.name)
    next!(validation_progress)

    # filter for examples w/ target data AND structure input
    filter!(r -> r.name ∈ found_cif, df) ##! SLOW
    next!(validation_progress)

    @assert df.name != [] "$target_csv and $xtal_dir have no structure names in common."
    next!(validation_progress)
    @info "$(length(df.name)) inputs."

    # sample list
    if args[:samples] ≠ 0
        @assert args[:samples] <= length(df.name) "Requested number of samples ($(args[:samples])) exceeds number of inputs."
        @info "Randomly sampling $(args[:samples])."
        df = df[in.(df.name, [sample(df.name, args[:samples], replace=false)]), :]
    end

    if args[:target] != ""
        @info "Checking target values."
        target_min = args[:target_lo]
        target_max = args[:target_hi]
        # filter xtal_list by 
        filter!(row -> row[args[:target]] >= target_min, df)
        filter!(row -> row[args[:target]] <= target_max, df)
        @assert length(df.name) > 0 "No target values within range [$target_min, $target_max]"
    end

    if args[:pvec]
        @everywhere args[:pvec_columns] = readlines(joinpath(rc[:paths][:data], "pvec_columns.txt"))
    end

    # make sure names in df are correct data type
    df.name = String.(df.name)

    # check that inputs can be read and start their xtal_data temp dicts
    loadcheck(df.name, args[:charge_tolerance], temp_dir)

    # process inputs to primitive cells
    xtals2primitive(df.name, temp_dir)

    # infer bonds, test quality
    good_xtals = bondNclassify(df.name, temp_dir)
    @assert length(good_xtals) > 0 "No structures passed bonding inspection."
    @info "$(length(good_xtals)) good bonded xtals."
    
    # determine atom node encoding scheme
    element_to_int, encoding_length = determine_encoding(good_xtals, args, temp_dir)
    @assert encoding_length > 0

    # final selection of target df rows
    filter!(row -> row.name in good_xtals, df)
    @assert length(df.name) == length(good_xtals)

    # process the graphs into ML inputs
    good = process_examples(good_xtals, element_to_int, df, args, temp_dir)

    good_xtals = good_xtals[good]

    dataset = consolidate_data(good_xtals, element_to_int, df, args, temp_dir)

    export_data(dataset, args)

    return good_xtals
end
