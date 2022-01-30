module Processing_test

using Test, MLMolGraph, JLD2, XtalsPyTools

@testset "xtals2primitive" begin
    clear_cache()
    xtal_name = "str_m2_o10_o29_pcu_sym.138.cif"
    xtal = Crystal(xtal_name)
    @test xtal.atoms.n == 92
    xtals2primitive([xtal_name])
    @load joinpath(rc[:cache][:primitive], xtal_name) obj
    xtal = obj
    @test xtal.atoms.n == 92
end

end