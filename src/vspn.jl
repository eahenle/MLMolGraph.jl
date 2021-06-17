struct VSPNConfig
    probe::Molecule
    ljff::LJForceField
end

VSPNConfig(probe::String, ljff::String) = VSPNConfig(Molecule(probe), LJForceField(ljff))


function show(config::VSPNConfig)
    print("VSPN Configuration\nForcefield: $(config.forcefield.name)\nProbe: $(config.molecule.species)")
end


function collect_atoms!(graph::MetaGraph, xtal::Crystal)
    @debug "Collecting atoms"  
    for i ∈ 1:xtal.atoms.n
        add_vertex!(graph, Dict(:type => :A, :species => xtal.atoms.species[i]))
    end
    for i ∈ 1:xtal.atoms.n
        for j ∈ neighbors(graph, i)
            add_edge!(graph, i, j, Dict(:distance => get_prop(xtal.bonds, i, j, :distance), :type => :AA))
        end
    end
end


function collect_vertices!(graph::MetaGraph, vt::VoroTess)
    @debug "Collecting Voronoi vertices"
    points = unique_voro_pts(vt)
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


function remove_high_E_verts!(graph::MetaGraph, vt::VoroTess, config::VSPNConfig)
    @debug "Removing high-energy V nodes"
    high_E_verts = []
    for i ∈ (vt.xtal.atoms.n + 1):nv(graph) # loop over vertex indices
        # translate probe to vertex
        translate_to!(config.probe,
            shift_back(get_prop(graph, i, point.coords), vt.xtal.box), # needs to be in the Xtals reference frame
            vt.xtal.box)
        # test the Van der Waals potential
        if vdw_energy(vt.xtal, config.probe, config.ljff) > 0
            append!(high_E_verts, i)
        end
    end
    # cut out the vertices w/ high potential energies
    rem_vertices!(graph, high_E_verts)
end


function mark_radii!(graph::MetaGraph, vt::VoroTess)
    @debug "Marking V node radii"
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


function greedy_selection1!(keep_vertices::Vector{Int}, radius_order::Vector{Int}, vt::VoroTess, graph::MetaGraph)
    @debug "First round of greedy algorithm"
    for i ∈ (vt.xtal.atoms.n .+ radius_order):nv(graph) # loop over vertices, largest radii first
        overlapping = false
        i_x = get_prop(graph, i, :point).coords
        i_r = get_prop(graph, i, :radius)
        for j ∈ findall(keep_vertices) # loop over vertices being kept and check if i is within radius of j ## TODO only check within local polytopes?
            if pbc_distance(i_x, get_prop(graph, radius_order[j], :point).coords, vt.xtal.box) < i_r + get_prop(graph, radius_order[j], :radius)
                overlapping = true # circumspheres of i and j overlap, so don't add i
                break
            end
        end
        if overlapping
            keep_vertices[i] = true
        end
    end
end


function greedy_selection2!(keep_vertices::Vector{Int}, radius_order::Vector{Int}, vt::VoroTess, graph::MetaGraph)
    @debug "Second round of greedy algorithm"
    for i ∈ (vt.xtal.atoms.n .+ radius_order):nv(graph) # loop over vertices, largest radii first
        contained = false
        i_x = get_prop(graph, i, :point).coords
        i_r = get_prop(graph, i, :radius)
        for j ∈ findall(keep_vertices) # loop over vertices being kept and check if i is within radius of j ## TODO only check within local polytopes?
            if pbc_distance(i_x, get_prop(graph, radius_order[j], :point).coords, vt.xtal.box) < i_r
                contained = true # i is inside j's sphere, so don't add i
                break
            end
        end
        if !contained
            keep_vertices[i] = true
        end
    end
end


function calculate_lenses!(graph::MetaGraph, vt::VoroTess)
    @debug "Calculating VV edge lenses"
    for i ∈ (vt.xtal.atoms.n + 1):nv(graph) # loop over remaining vertices
        i_r = get_prop(graph, i, :radius)
        i_x = get_prop(graph, i, :point).coords
        for j ∈ i:nv(graph) # loop over vertex pairs
            # make a ridge edge if radii overlap ## TODO only check within local polytopes?
            j_r = get_prop(graph, j, :radius)
            dist = pbc_distance(i_x, get_prop(graph, j, :point).coords, vt.xtal.box)
            if dist < j_r + i_r
                add_edge!(graph, i, j, Dict(:type => :VV, :lens => lens(i_r, i_j, dist)))
            end
        end
    end
end


function vspn_graph(xtal::Crystal, config::VSPNConfig)::MetaGraph
    @debug "Getting VSPN graph"
    # get voronoi tesselation
    vt = voronoi_tesselation(xtal)
    @debug "Got tesselation"
    # merge graphs
    output_graph = MetaGraph()
    # atoms and bonds
    collect_atoms!(output_graph, xtal)
    @debug "Collected atoms"
    # Voronoi vertices
    collect_vertices!(output_graph, vt)
    @debug "Collected vertices" 
    # remove vertices (pass 1)
    remove_high_E_verts!(output_graph, vt, config)
    @debug "Removed high-energy vertices" 
    # find sphere radii
    mark_radii!(output_graph, vt)
    # sort vertices by radius, descending order
    radius_order = sortperm([get_prop(output_graph, i, :radius) for i ∈ (vt.nb_atoms + 1):nv(output_graph)], rev=true)
    @debug "Marked and sorted V node radii"  
    # determine which vertices to keep
    keep_vertices = falses(length(radius_order))
    # tag vertices for keeping (pass 2)
    greedy_selection1!(keep_vertices, radius_order, vt, output_graph)
    @debug "Did first greedy selection" 
    # tag vertices for keeping (pass 3)
    greedy_selection2!(keep_vertices, radius_order, vt, output_graph)
    @debug "Did second greedy selection" 
    # remove vertices (passes 2 & 3)
    rem_vertices!(output_graph, ((vt.nb_atoms + 1):nv(output_graph))[.!keep_vertices])
    @debug "Filtered V nodes" 
    # make VV edges and label by lens volume
    calculate_lenses!(output_graph, vt)
    @debug "Calculated lenses" 
    return output_graph
end

vspn_graph(xtal::String, probe::String, ljff::String) = vspn_graph(Crystal(xtal), VSPNConfig(probe, ljff))
