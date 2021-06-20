module Run_process_test

using Test, MLMolGraph

# Note: must specify all args (argparse not applied)

@testset "vspn" begin
    clear_cache()
    run_process(Dict(
        :target => Symbol("deliverable capacity [v STP/v]"),
        :data => joinpath(pwd(), "data"),
        :cache => joinpath(pwd(), "data", "cache"),
        :bonds => false,
        :angles => false,
        :vspn => true,
        :probe => "CH4",
        :forcefield => "UFF"
    ))
    @test true
end

@testset "vspn + angles" begin
    clear_cache()
    run_process(Dict(
        :target => Symbol("deliverable capacity [v STP/v]"),
        :data => joinpath(pwd(), "data"),
        :cache => joinpath(pwd(), "data", "cache"),
        :bonds => false,
        :angles => true,
        :vspn => true,
        :probe => "CH4",
        :forcefield => "UFF"
    ))
    @test true
end

@testset "bonds" begin
    clear_cache()
    run_process(Dict(
        :target => Symbol("deliverable capacity [v STP/v]"),
        :data => joinpath(pwd(), "data"),
        :cache => joinpath(pwd(), "data", "cache"),
        :bonds => true,
        :angles => false,
        :vspn => false,
        :probe => "",
        :forcefield => ""
    ))
    @test true
end

@testset "angles" begin
    clear_cache()
    run_process(Dict(
        :target => Symbol("deliverable capacity [v STP/v]"),
        :data => joinpath(pwd(), "data"),
        :cache => joinpath(pwd(), "data", "cache"),
        :bonds => true,
        :angles => true,
        :vspn => false,
        :probe => "",
        :forcefield => ""
    ))
    @test true
end

end
