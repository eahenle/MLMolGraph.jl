struct VoroTess
    xtal
    cells
    box
    atom_pts
    n_dict
end


struct VoroPoint
    coords
    cells
end


struct VoroCell
    vertices
    neighbors
end


function neighbor_dict(n_list)
    n_dict = Dict()
    for (i, j) ∈ n_list
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
function voronoi_tesselation(points::Matrix{Float64}, box_params::Array{Float64,1})
    freud = rc[:freud]
    box = freud.box.Box.from_box(box_params)
    voro = freud.locality.Voronoi()
    voro.compute((box, points))
    cells = voro.polytopes
    @info "FOO" voro voro.nlist voro.nlist.point_indices
    n_dict = neighbor_dict(voro.nlist.point_indices)
    return VoroTess(xtal, cells, box, points, n_dict)
end

voronoi_tesselation(xtal::Crystal) = voronoi_tesselation(shift_coords(xtal), convert_box(xtal.box))


function unique_voro_pts(vt)
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