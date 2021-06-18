module Processing_test

using Test, MLMolGraph

@testset "xtals2primitive" begin
    xtal_name = "linker91_CH2_linker11_NH_fxt_relaxed.cif"
    xtal = Crystal(xtal_name)
    @test xtal.atoms.n = 1728
    xtals2primitive([xtal_name])
    old_xtal_path = rc[:paths][:crystals]
    rc[:paths][:crystals] = rc[:cache][:primitive]
    xtal = Crystal(xtal_name)
    @test xtal.atoms.n = 576
    rc[:paths][:crystals] = old_xtal_path
end

end