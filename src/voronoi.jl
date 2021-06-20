struct VoroTess
    xtal::Crystal
    cells::Vector{Matrix{Float64}}
    box_params::Vector{Float64}
    atom_pts::Matrix{Float64}
    n_dict::Dict{Int,Vector{Int}}
end


struct VoroPoint
    coords::Vector{Float64}
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
    print("Crystal: $(vt.xtal.name)")
    print("$(length(vt.cells)) cells")
    print("Box parameters: $(vt.box_params)")
    print("$(size(vt.atom_pts, 2)) atom points")
    print("Neighbors: $(vt.n_dict)")
end
display(vt::VoroTess) = show(vt)


function neighbor_dict(A::Vector{Int}, B::Vector{Int})::Dict{Int,Vector{Int}}
    @assert length(A) == length(B)
    n_dict = Dict()
    for k ∈ 1:length(A)
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
    # round point coords to 6 digits for uniqueness testing
    for i ∈ 1:length(cells)
        for j ∈ 1:length(cells[i])
            cells[i][j] = round(cells[i][j], digits=6)
        end
    end
    # correct for Python indexing
    A = [i+1 for i ∈ voro.nlist.query_point_indices]
    B = [i+1 for i ∈ voro.nlist.point_indices]
    n_dict = neighbor_dict(A, B)
    return VoroTess(xtal, cells, box_params, points, n_dict)
end

voronoi_tesselation(xtal::Crystal) = voronoi_tesselation(xtal, shift_coords(xtal), convert_box(xtal.box))


function unique_voro_pts(vt::VoroTess)::Vector{VoroPoint}
    points = VoroPoint[]
    cells = Dict{Int,VoroCell}()
    # loop over vt.cells
    for (c, cell) ∈ enumerate(vt.cells)
        # some cells have no entry in the neighbor list of the tesselation; skip them
        if c ∉ keys(vt.n_dict)
            continue
        end
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
        cells[c] = VoroCell(cell[unique_col_ids,:], vt.n_dict[c])
    end
    # build a DataFrame w/ each row a single vertex of a single cell (merge them later)
    df = DataFrame(:cell_index => Int[], :pt_coords => Vector{Float64}[])
    # loop over cells
    for (c, cell) ∈ cells
        # loop over points
        for i ∈ 1:size(cell.vertices,1)
            point = cell.vertices[i,:]
            # add point coords, cell index to row of df
            append!(df, DataFrame(:cell_index => [c], :pt_coords => [point]))
        end
    end
    # group df by coords
    gdf = groupby(df, :pt_coords)
    # combine grouped rows by merging cell indices into list
    for group ∈ gdf
        append!(points, [VoroPoint(group.pt_coords[1], unique(group.cell_index))])
    end
    return points
end