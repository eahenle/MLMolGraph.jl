module Voronoi_test

using Test, MLMolGraph

@testset "unique_voro_pts" begin
    xtal = Crystal("linker91_CH2_linker1_NH_hcb_relaxed.cif", remove_duplicates = true)
    vt = voronoi_tesselation(primitive_cell(xtal))
    unique_points = MLMolGraph.unique_voro_pts(vt)
    @test true
end

end