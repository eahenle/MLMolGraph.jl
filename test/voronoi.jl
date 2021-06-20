module Voronoi_test

using Test, MLMolGraph

@testset "n_dict" begin
    A = [1, 1, 2, 2, 3, 3, 4, 4]
    B = [2, 3, 1, 4, 1, 4, 2, 3]
    n_dict = MLMolGraph.neighbor_dict(A, B)
    @test n_dict[1] == [2,3] 

    xtal = Crystal(
        "square grid",
        unit_cube(),
        Atoms(
            [:H, :H, :H, :H],
            Frac([
                0.25 0.75 0.25 0.75;
                0.75 0.75 0.25 0.25;
                0.00 0.00 0.00 0.00
            ])
        ),
        Charges{Frac}(0)
    )
    vt = voronoi_tesselation(xtal)
    @test vt.n_dict[1] == [1,2,3] # ONLY RECOGNIZES FACET-SHARING, SO NO 1 → 4
end

@testset "unique_voro_pts" begin
    xtal = Crystal("linker91_CH2_linker1_NH_hcb_relaxed.cif", remove_duplicates = true)
    vt = voronoi_tesselation(primitive_cell(xtal))
    unique_points = MLMolGraph.unique_voro_pts(vt)
    pass = true
    for (i, vpi) ∈ enumerate(unique_points)
        for (j, vpj) ∈ enumerate(unique_points)
            if i == j
                continue
            end
            if isapprox(vpi.coords, vpj.coords)
                pass = false
                @error "Test failure" i j vpi.coords vpj.coords
                break
            end
        end
    end
    @test pass
end

end