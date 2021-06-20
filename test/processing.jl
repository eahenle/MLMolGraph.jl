module Processing_test

using Test, MLMolGraph, JLD2

@testset "xtals2primitive" begin
    xtal_name = "linker91_CH2_linker11_NH_fxt_relaxed.cif"
    xtal = Crystal(xtal_name)
    @test xtal.atoms.n == 1728
    xtals2primitive([xtal_name])
    @load joinpath(rc[:cache][:primitive], xtal_name) obj
    xtal = obj
    @test xtal.atoms.n == 576
end

end