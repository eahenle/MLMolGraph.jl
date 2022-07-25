module Processing_test

using Graphs, Test, MLMolGraph, XtalsPyTools

temp_dir = joinpath(pwd(), ".temp")
if !isdir(temp_dir)
    mkpath(temp_dir)
end

@testset "xtals2primitive" begin
    xtal_name = "str_m3_o2_o18_pcu_sym.41"
    MLMolGraph.loadcheck([xtal_name], 1e-4, temp_dir)
    MLMolGraph.xtals2primitive([xtal_name], temp_dir)
    xtal = MLMolGraph.load_data(xtal_name, temp_dir)[:primitive_cell]
    @test xtal.atoms.n == 90
end

end