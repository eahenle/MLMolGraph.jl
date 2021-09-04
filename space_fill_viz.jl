import Pkg
Pkg.activate(".")

using JLD2, LightGraphs, MetaGraphs, MLMolGraph
import MLMolGraph.VoroPoint


function write_vis_data(xtal::Crystal, voro_graph::MetaGraph, label::String; color="He")
    # extract Voronoi vertices
    v_vert_idx = (vt.xtal.atoms.n + 1):nv(voro_graph)
    points = [get_prop(voro_graph, i, :point).coords for i in v_vert_idx]

    # Cartesian coordinates and radii
    r = [get_prop(voro_graph, i, :radius) for i in v_vert_idx]
    x = [p[1] for p in points]
    y = [p[2] for p in points]
    z = [p[3] for p in points]
    
    # format into XYZ format for VisIt
    spheres = "$(length(points))\n#\n"
    for i in 1:length(points)
        spheres *= "$color\t$(x[i])\t$(y[i])\t$(z[i])\t$(r[i])\n"
    end

    # write outputs
    if !isdir("spacefill")
        mkdir("spacefill")
    end
    write_xyz(xtal, "spacefill/atoms.xyz")
    write_vtk(xtal.box, "spacefill/cell.vtk")
    write_bond_information(xtal, "spacefill/bonds.vtk")
    ## TODO add V-V edges
    open("spacefill/$label-spheres.xyz", "w") do io
        print(io, spheres)
    end
end


@info "Loading crystal"
# load example crystal and get bonds
xtal = load_object("data/cache/primitive/str_m3_o23_o24_pcu_sym.19.cif")
infer_bonds!(xtal, true)

@info "Finding Voronoi tesselation"
# perform VT on atoms
vt = voronoi_tesselation(xtal)

@info "Building graph"
# build the Voro graph and label the sphere radii
voro_graph = MetaGraph()
MLMolGraph.collect_atoms!(voro_graph, xtal)
MLMolGraph.collect_vertices!(voro_graph, vt)
@info "Finding accessible space"
MLMolGraph.remove_high_E_verts!(voro_graph, vt, MLMolGraph.VSPNConfig("CH4", "UFF"))
MLMolGraph.mark_radii!(voro_graph, vt)

@info "Sorting by included sphere radius"
# sort V nodes by radius (descending)
radius_order = sortperm([get_prop(voro_graph, i, :radius) for i ∈ (vt.xtal.atoms.n + 1):nv(voro_graph)], rev=true)
keep_vertices = falses(length(radius_order))

@info "First pass of space-filling algorithm"
# take first pass of greedy algorithm, record "good" nodes
MLMolGraph.greedy_selection1!(keep_vertices, radius_order, vt, voro_graph)
pass1_keeps = deepcopy(keep_vertices)
pass1_graph = deepcopy(voro_graph)
drop_verts = sort(((vt.xtal.atoms.n + 1):nv(pass1_graph))[.! pass1_keeps], rev=true)
for v ∈ drop_verts
    rem_vertex!(pass1_graph, v)
end
write_vis_data(xtal, pass1_graph, "pass1")

@info "Second pass of space-filling algorithm"
# take second pass of greedy algorithm, record "good" nodes
MLMolGraph.greedy_selection2!(keep_vertices, radius_order, vt, voro_graph)
pass2_keeps = deepcopy(keep_vertices)
pass2_graph = deepcopy(voro_graph)
drop_verts = sort(((vt.xtal.atoms.n + 1):nv(pass2_graph))[.! pass2_keeps], rev=true)
for v ∈ drop_verts
    rem_vertex!(pass2_graph, v)
end
write_vis_data(xtal, pass2_graph, "pass2")

@info "Calculating lenses"
# calculate the V-V edge weights
MLMolGraph.calculate_lenses!(pass2_graph, vt)
## TODO write lens data

@info "Done."
