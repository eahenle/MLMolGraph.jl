using FIGlet, Test, MLMolGraph

FIGlet.render("MLMolGraph", FIGlet.availablefonts()[64])

@testset "bonds" begin
    clear_cache()
    run_process(Dict(
        :target => "deliverable capacity [v STP/v]",
        :data => joinpath(pwd(), "data"),
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
        :target => "deliverable capacity [v STP/v]",
        :data => joinpath(pwd(), "data"),
        :bonds => true,
        :angles => true,
        :vspn => false,
        :probe => "",
        :forcefield => ""
    ))
    @test true
end

@testset "vspn" begin
    clear_cache()
    run_process(Dict(
        :target => "deliverable capacity [v STP/v]",
        :data => joinpath(pwd(), "data"),
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
        :target => "deliverable capacity [v STP/v]",
        :data => joinpath(pwd(), "data"),
        :bonds => false,
        :angles => true,
        :vspn => true,
        :probe => "CH4",
        :forcefield => "UFF"
    ))
    @test true
end
