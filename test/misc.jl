module Misc_test

using Test, MLMolGraph

@testset "PBC distance" begin
    points = [ # grid of 0.5x0.5x0.5 cube vertices
        0.25 0.25 0.25 0.75 0.75 0.75 0.25 0.75;
        0.25 0.25 0.75 0.25 0.25 0.75 0.75 0.75;
        0.25 0.75 0.25 0.25 0.75 0.25 0.75 0.75
    ]
    cubic_box = unit_cube()
    @test MLMolGraph.pbc_distance(points[:,1], points[:,2], cubic_box) == 0.5
    @test MLMolGraph.pbc_distance(points[:,1], points[:,8], cubic_box) == √(0.75)

    triclinic_box = Box(1., 1., 1., π/3, π/3, π/3)
    pass = true
    for i ∈ 1:8
        for j ∈ 1:8
            pt_i = points[:,i]
            pt_j = points[:,j]
            cb_dist = MLMolGraph.pbc_distance(pt_i, pt_j, cubic_box)
            tc_dist = MLMolGraph.pbc_distance(pt_i, pt_j, triclinic_box)
            if cb_dist ≠ tc_dist
                @error "Test failure" i j pt_i pt_j
            end
        end
    end
    @test pass
end

end
