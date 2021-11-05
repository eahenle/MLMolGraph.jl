struct VSPN_Input_Struct
    xtal_name::String
    adj_mat::SparseArrays.SparseMatrixCSC
    X::SparseArrays.SparseMatrixCSC
    y::Float64
end


struct VSPNConfig
    probe::Molecule
    ljff::LJForceField
end

VSPNConfig(probe::String, ljff::String) = VSPNConfig(Molecule(probe), LJForceField(ljff))


function show(config::VSPNConfig)
    print("VSPN Configuration\nForcefield: $(config.forcefield.name)\nProbe: $(config.molecule.species)")
end


function collect_atoms!(graph::MetaGraph, xtal::Crystal)
    for i ∈ 1:xtal.atoms.n
        add_vertex!(graph)
        set_props!(graph, i, Dict(:type => :A, :species => xtal.atoms.species[i]))
    end
    for i ∈ 1:xtal.atoms.n
        for j ∈ neighbors(graph, i)
            add_edge!(graph, i, j)
            set_props!(graph, i, j, Dict(:distance => get_prop(xtal.bonds, i, j, :distance), :type => :AA))
        end
    end
end


function collect_vertices!(graph::MetaGraph, vt::VoroTess, args)
    points = unique_voro_pts(vt)
    if args[:verbose]
        @info "$(length(points)) unique voro pts"
    end
    for point ∈ points
        add_vertex!(graph, Dict(:type => :V, :point => point))
    end
    # make bridge edges to neighbor atoms, annotating distance
    for i ∈ (vt.xtal.atoms.n + 1):(nv(graph)) # loop over vertex indices
        point = get_prop(graph, i, :point)
        for n ∈ point.cells # loop over indices (atoms) of adjacent cells
            add_edge!(graph, i, n, Dict(:type => :AV, :distance => pbc_distance(vt.atom_pts[n, :], point.coords, vt.box_params)))
        end
    end
end


function remove_high_E_verts!(graph::MetaGraph, vt::VoroTess, config::VSPNConfig, args)
    high_E_verts = Int[]
    for i ∈ (vt.xtal.atoms.n + 1):nv(graph) # loop over vertex indices
        # translate probe to vertex
        translate_to!(config.probe, # needs to be in the Xtals reference frame
            Cart(shift_back(get_prop(graph, i, :point).coords, vt.xtal.box)))
        # test the Van der Waals potential
        if vdw_energy(vt.xtal, config.probe, config.ljff) > 0
            append!(high_E_verts, i)
        end
    end
    # cut out the vertices w/ high potential energies
    for i ∈ sort(high_E_verts, rev=true) ## * BEWARE this BREAKS CORRESPONDENCE in arrays outside of graph!!
        rem_vertex!(graph, i)
    end

    if args[:verbose]
        n = count(v -> get_prop(graph, v, :type) == :V, 1:nv(graph))
        @info "$n type-V nodes after high-E removal"
    end
end


function mark_radii!(graph::MetaGraph, vt::VoroTess, args)
    for i ∈ (vt.xtal.atoms.n + 1):nv(graph) # loop over vertex indices
        min_dist = Inf
        for n ∈ neighbors(graph, i) # loop over adjacent cell atoms
            if get_prop(graph, n, :type) == :A
                dist = get_prop(graph, i, n, :distance)
                if dist < min_dist
                    min_dist = dist # find shortest bridge edge
                end 
            end
        end
        set_prop!(graph, i, :radius, min_dist)
    end
end


function greedy_selection1!(keep_vertices::BitVector, radius_order::Vector{Int}, vt::VoroTess, voro_graph::MetaGraph)
    for i ∈ radius_order # loop over vertices, largest radii first
        overlapping = false
        i_x = get_prop(voro_graph, i, :point).coords
        i_r = get_prop(voro_graph, i, :radius)
        for j ∈ findall(keep_vertices) # loop over vertices being kept and check if i is within radius of j ## TODO only check within local polytopes?
            j_r = get_prop(voro_graph, j, :radius)
            if pbc_distance(i_x, get_prop(voro_graph, j, :point).coords, vt.xtal.box) < i_r + j_r
                overlapping = true # circumspheres of i and j overlap and are dissimilar, so don't add i
                break
            end
        end
        if !overlapping
            keep_vertices[i] = true
        end
    end
end


function greedy_selection2!(keep_vertices::BitVector, radius_order::Vector{Int}, vt::VoroTess, voro_graph::MetaGraph)
    for i ∈ radius_order # loop over vertices, largest radii first
        contained = false
        i_x = get_prop(voro_graph, i, :point).coords
        for j ∈ findall(keep_vertices) # loop over vertices being kept and check if i is within radius of j ## TODO only check within local polytopes?
            d = pbc_distance(i_x, get_prop(voro_graph, radius_order[j], :point).coords, vt.xtal.box)
            if d < minimum([get_prop(voro_graph, i, :radius), get_prop(voro_graph, j, :radius)])
                contained = true # i is inside j's sphere, so don't add i
                break
            end
        end
        if !contained
            keep_vertices[i] = true
        end
    end
end


function calculate_lenses!(voro_graph::MetaGraph, vt::VoroTess)
    for i ∈ 1:nv(voro_graph) # loop over remaining vertices
        i_r = get_prop(voro_graph, i, :radius)
        i_x = get_prop(voro_graph, i, :point).coords
        for j ∈ i:nv(voro_graph) # loop over vertex pairs
            # make a ridge edge if radii overlap ## TODO only check within local polytopes?
            j_r = get_prop(voro_graph, j, :radius)
            dist = pbc_distance(i_x, get_prop(voro_graph, j, :point).coords, vt.xtal.box)
            if dist < j_r + i_r
                add_edge!(voro_graph, i, j, Dict(:type => :VV, :lens => lens(i_r, j_r, dist)))
            end
        end
    end
end


function extract_voro_graph(combined_graph::MetaGraph, vt::VoroTess)::MetaGraph
    voro_graph = MetaGraph()
    # pull out the V-nodes with all their meta-data
    for v in (vt.xtal.atoms.n+1):nv(combined_graph)
        add_vertex!(voro_graph, props(combined_graph, v))
    end
    # re-create edges
    for e in edges(combined_graph)
        add_edge!(voro_graph, src(e)-vt.xtal.atoms.n, dst(e)-vt.xtal.atoms.n, props(combined_graph, e))
    end
    return voro_graph
end


function vspn_graph(xtal::Crystal, config::VSPNConfig, args)::MetaGraph
    # get voronoi tesselation
    vt = voronoi_tesselation(xtal)
    # merge graphs
    combined_graph = MetaGraph()
    # atoms and bonds
    collect_atoms!(combined_graph, xtal)
    # Voronoi vertices
    collect_vertices!(combined_graph, vt, args)
    # remove vertices (pass 1)
    remove_high_E_verts!(combined_graph, vt, config, args)
    # find sphere radii
    mark_radii!(combined_graph, vt, args)
    # separate the bonding and Voronoi graphs
    voro_graph = extract_voro_graph(combined_graph, vt) ## ! only using voro_graph; will need refactor to allow combined_graph
    # sort vertices by radius, descending order
    radius_order = sortperm([get_prop(voro_graph, i, :radius) for i ∈ 1:nv(voro_graph)], rev=true)
    # determine which vertices to keep
    keep_vertices = falses(length(radius_order))
    # tag vertices for keeping (pass 2)
    greedy_selection1!(keep_vertices, radius_order, vt, voro_graph)
    if args[:verbose]
        @info "$(sum(keep_vertices)) type-V nodes selected after greedy round 1"
    end
    # tag vertices for keeping (pass 3)
    greedy_selection2!(keep_vertices, radius_order, vt, voro_graph)
    if args[:verbose]
        @info "$(sum(keep_vertices)) type-V nodes selected after greedy round 2"
    end
    # remove vertices (passes 2 & 3)
    drop_verts = sort((1:nv(voro_graph))[.!keep_vertices], rev=true)
    for v ∈ drop_verts
        rem_vertex!(voro_graph, v)
    end
    # make VV edges and label by lens volume
    calculate_lenses!(voro_graph, vt)
    return voro_graph
end

vspn_graph(xtal::String, probe::String, ljff::String, args) = vspn_graph(Crystal(xtal), VSPNConfig(probe, ljff), args)
