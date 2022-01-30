module Run_process_test

using Test, MLMolGraph, NPZ

# Note: must specify all args (argparse not applied)
function make_arg_dict(target::String)::Dict{Symbol, Any}
    Dict(
        :target => Symbol(target),
        :data => joinpath(pwd(), "data"),
        :cache => joinpath(pwd(), "data", "cache"),
        :bonds => false,
        :angles => false,
        :vspn => false,
        :probe => "",
        :forcefield => "",
        :samples => 0,
        :primitive => true,
        :tolerance => 1e-5,
        :verbose => true,
        :one_hot => "n",
        :vectors => false
    )
end

@testset "bonds" begin
    clear_cache()
    arg_dict = make_arg_dict("working_capacity_vacuum_swing [mmol/g]")
    arg_dict[:bonds] = true
    xtal_names = run_process(arg_dict)
    for xtal_name in xtal_names
        xtal_name = split(xtal_name, ".cif")[1]
        X = npzread(joinpath(rc[:paths][:graphs], xtal_name * "_node_features.npy"))
        # test that standard encoding rows all sum to 1
        @assert all(1 .== [sum(X[i, :]) for i in 1:size(X, 1)])
    end
    @test length(xtal_names) > 0
end

@testset "angles" begin
    clear_cache()
    arg_dict = make_arg_dict("working_capacity_vacuum_swing [mmol/g]")
    arg_dict[:angles] = true
    xtal_names = run_process(arg_dict)
    @test length(xtal_names) > 0
end

@testset "vspn" begin
    clear_cache()
    arg_dict = make_arg_dict("working_capacity_vacuum_swing [mmol/g]")
    arg_dict[:vspn] = true
    arg_dict[:probe] = "CH4"
    arg_dict[:forcefield] = "UFF"
    xtal_names = run_process(arg_dict)
    @test length(xtal_names) > 0
end

end
