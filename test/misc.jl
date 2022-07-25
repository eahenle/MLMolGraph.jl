module Misc_test

using Test, MLMolGraph
import MLMolGraph.pbc_distance

@testset "PBC distance" begin
    points = [ # grid of 0.5x0.5x0.5 cube vertices
        0.25 0.25 0.25 0.75 0.75 0.75 0.25 0.75;
        0.25 0.25 0.75 0.25 0.25 0.75 0.75 0.75;
        0.25 0.75 0.25 0.25 0.75 0.25 0.75 0.75
    ]
    cubic_box = unit_cube()
    @test pbc_distance(points[:,1], points[:,2], cubic_box) == 0.5
    @test pbc_distance(points[:,1], points[:,8], cubic_box) == √(0.75)

    triclinic_box = Box(1., 1., 1., π/3, π/3, π/3)
    @test all([pbc_distance(pt_i, pt_j, cubic_box) == pbc_distance(pt_j, pt_i, triclinic_box) for i in 1:8 for j in 1:8])
end

end
