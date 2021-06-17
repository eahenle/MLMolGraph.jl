struct VoroTess
    xtal::Crystal
    cells::Vector{Matrix{Float64}}
    box_params::Vector{Float64}
    atom_pts::Matrix{Float64}
    n_dict::Dict{Int,Vector{Int}}
end


struct VoroPoint
    coords::Matrix{Float64}
    cells::Vector{Int}
end


struct VoroCell
    vertices::Matrix{Float64}
    neighbors::Vector{Int}
end


function show(vc::VoroCell)
    print("Vertices: $(vc.vertices)")
    print("Neighboring cells: $(vc.neighbors)")
end


function show(vp::VoroPoint)
    print("Coords: $(vp.coords)")
    print("Vertex of cells: $(vp.cells)")
end


function show(vt::VoroTess)
    print("Xtal: $(vt.xtal.name)")
    print("$(length(vt.cells)) cells")
    print("Box parameters: $(vt.box_params)")
    print("$(size(vt.atom_pts, 2)) atom points")
    print("Neighbors: $(vt.n_dict)")
end
display(vt::VoroTess) = show(vt)


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
    # correct for Python indexing
    A = [i+1 for i ∈ voro.nlist.query_point_indices]
    B = [i+1 for i ∈ voro.nlist.point_indices]
    n = voro.nlist.num_points
    n_dict = neighbor_dict(A, B, n)
    return VoroTess(xtal, cells, box_params, points, n_dict)
end

voronoi_tesselation(xtal::Crystal) = voronoi_tesselation(xtal, shift_coords(xtal), convert_box(xtal.box))


function unique_voro_pts(vt::VoroTess)::Vector{VoroPoint}
    @debug "Finding unique polytope vertices" vt.cells
    points = VoroPoint[]
    cells = Array{VoroCell}(undef, length(vt.cells))
    # loop over vt.cells
    for (c, cell) ∈ enumerate(vt.cells)
        # each entry in vt.cells is a matrix of point coordinates, but transposed from what's best
        mat = cell'
        # add the unique columns from the matrix to the cell struct
        unique_col_ids = Int[]
        unique_cols = Vector{Float64}[]
        for i ∈ 1:size(mat,2) # loop over columns
            col = mat[:,i]
            if col ∉ unique_cols
                append!(unique_cols, [col])
                append!(unique_col_ids, i)
            end
        end
        # pull the list of neighboring cell indices and make the cell struct
        @debug "Making VoroCell" c 
        @debug vt.n_dict
        @debug cell[unique_col_ids,:]
        cells[c] = VoroCell(cell[unique_col_ids,:], vt.n_dict[c])
    end
    df = DataFrame(:cell_index => Int[], :pt_coords => Matrix{Float64}[])
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