# All args must be set for tests to pass. Keep this Dict updated w/ new args from table.
DEFAULT_ARGS = Dict(
    # `store_true` args (must be false by default)
    :vectors        => false,
    :bonds          => false,
    :angles         => false,
    :verbose        => false,
    :pvec           => false,

    # string args
    :data               => rc[:paths][:data],
    :target             => "",
    :dataset_output     => "_DATASET_.PKL",
    
    # int args
    :seed               => Int(floor(rand() * 1e3)),
    :samples            => 0,
    
    # float args
    :charge_tolerance   => 1e-5,
    :target_hi          => Inf,
    :target_lo          => -Inf
)

@add_arg_table! argparser begin
    "--angles"
        help = "write ML inputs for bonding angles"
        action = :store_true

    "--bonds"
        help = "write ML inputs for the bonding graph"
        action = :store_true

    "--seed"
        help = "set the random seed for reproducible random dataset sampling"
        arg_type = Int
        default = DEFAULT_ARGS[:seed]

    "--data"
        help = "root of data tree. inputs loaded from `crystals` subdirectory."
        arg_type = String
        default = DEFAULT_ARGS[:data]
    
    "--verbose"
        help = "whether or not to print extra @info output"
        action = :store_true

    "--vectors"
        help = "write bond vectors"
        action = :store_true

    "--charge_tolerance"
        help = "set net charge tolerance for Crystal()"
        arg_type = Float64
        default = DEFAULT_ARGS[:charge_tolerance]
         
    "--samples"
        help = "number of inputs to sample (0 -> use all)"
        arg_type = Int
        default = DEFAULT_ARGS[:samples]

    "--target"
        help = """name of target column for constraining values. if left blank, all targets will be included in the dataset["y"] tensor."""
        arg_type = String
        default = DEFAULT_ARGS[:target]

    "--target_hi"
        help = "upper bound for target value"
        arg_type = Float64
        default = DEFAULT_ARGS[:target_hi]

    "--target_lo"
        help = "lower bound for target value"
        arg_type = Float64
        default = DEFAULT_ARGS[:target_lo]

    "--dataset_output"
        help = "path to file where processed data will be saved"
        default = DEFAULT_ARGS[:dataset_output]

    "--pvec"
        help = "whether or not to look for physical vector information"
        default = DEFAULT_ARGS[:pvec]
end
