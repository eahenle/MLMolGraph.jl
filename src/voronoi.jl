struct VoroTess
    xtal::Crystal
    cells::Vector{Matrix{Float64}}
    box_params::Vector{Float64}
    atom_pts::Matrix{Float64}
    n_dict::Dict{Int,Int}
end


struct VoroPoint
    coords::Matrix{Float64}
    cells::Vector{Int}
end


struct VoroCell
    vertices::Vector{Int}
    neighbors::Vector{Int}
end


function neighbor_dict(A::Vector{Int}, B::Vector{Int}, n::Int)::Dict{Int,Vector{Int}}
    n_dict = Dict()
    for k ∈ 1:n
        i = A[k]
        j = B[k]
        if i ∈ keys(n_dict)
            append!(n_dict[i], j)
        else
            n_dict[i] = [j]
        end
        if j ∈ keys(n_dict)
            append!(n_dict[j], i)
        else
            n_dict[j] = [i]
        end
    end
    for (key, value) ∈ n_dict
        unique!(n_dict[key])
    end
    return n_dict
end


# calculates the Voronoi tesselation
function voronoi_tesselation(xtal::Crystal, points::Matrix{Float64}, box_params::Array{Float64,1})::VoroTess
    freud = rc[:freud]
    box = freud.box.Box.from_box(box_params)
    voro = freud.locality.Voronoi()
    voro.compute((box, points))
    cells = voro.polytopes
    A = [i for i ∈ voro.nlist.query_point_indices]
    B = [i for i ∈ voro.nlist.point_indices]
    n = voro.nlist.num_points
    n_dict = neighbor_dict(A, B, n)
    return VoroTess(xtal, cells, box, points, n_dict)
end

voronoi_tesselation(xtal::Crystal) = voronoi_tesselation(xtal, shift_coords(xtal), convert_box(xtal.box))


function unique_voro_pts(vt::VoroTess)::Vector{VoroPoint}
    points = VoroPoint[]
    cells = Array{VoroCell}(undef, length(vt.cells))
    # loop over vt.cells
    for (c, cell) ∈ enumerate(vt.cells)
        # each entry in vt.cells is a list of points
        # add the unique points from the list to the cell struct
        # add the list of neighboring cells to the cell struct
        cells[c] = VoroCell(unique(cell), vt.n_dict[c])
    end
    df = DataFrame(:cell_index => [], :pt_coords => [])
    # loop over cells
    for (c, cell) ∈ enumerate(cells)
        # loop over points
        for point ∈ cell.vertices
            # add point coords, cell index to row of df
            append!(df, DataFrame(:cell_index => c, :pt_coords => point))
        end
    end
    # group df by coords
    gdf = groupby(df, :pt_coords)
    # combine grouped rows by merging cell indices into list
    cdf = combine(gdf, :cell_index => unique)
    # translate combined rows into VoroPoints
    for row ∈ rows(cdf)
        append!(points, VoroPoint(row.pt_coords, row.cell_index_unique))
    end
    return points
end