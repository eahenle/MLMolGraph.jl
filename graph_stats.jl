### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 9ee1c6ad-8660-422c-a1f5-d2b6b93875b3
begin # dev-hack
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 1e19d039-b180-409b-af52-610bec0fc169
using Graphs, JLD2, MetaGraphs, MLMolGraph, PorousMaterials, StatsBase

# ╔═╡ ea41da3c-53a4-4d60-8aa1-606903c3c3ae
begin
	voro_dir = "data/cache/vspn"
	bond_dir = "data/cache/bonded_xtals"
	nb_graphs = 100
end;

# ╔═╡ 6f89aa61-cb1c-4c55-b01f-1099e9a59376
xtal_names = String.(reduce(vcat, 
	split.(sample(readdir(voro_dir), nb_graphs, replace=false), ".jld2",
		keepempty=false)
));

# ╔═╡ 04970855-74ce-47fc-9642-0edf4563c9ca
function load_voro_graph(xtal_name::String; dir::String="")
    @load joinpath(dir, "$xtal_name.jld2") obj
    voro_graph, _ = obj
    return voro_graph
end;

# ╔═╡ 0239bc4c-741c-4176-a66c-526398d0fd44
voro_graphs = load_voro_graph.(xtal_names, dir=voro_dir);

# ╔═╡ f1f0f0f9-8e8d-445b-aa53-cb9ba1354ad3
function load_xtal(xtal_name::String; dir::String="")
	@load joinpath(dir, "$xtal_name.cif") obj
	xtal = obj[1]
	return xtal
end;

# ╔═╡ 8b1d1c34-46fc-4e75-8a1a-aeb3f19be118
crystals = load_xtal.(xtal_names, dir=bond_dir);

# ╔═╡ cd76fb5e-0a52-4655-8ff2-a1ec580cb68d
degree_population = reduce(vcat, degree.(voro_graphs));

# ╔═╡ da66a716-cf5e-4d51-ace2-ea3c671cefa0
begin
	using PyPlot
	figure()
	title("Voro-Graph Degree Distribution")
	xlabel("Node Degree")
	ylabel("Count")
	hist(degree_population, bins=maximum(degree_population)+1)
	gcf()
end

# ╔═╡ 028398e0-30a4-4ff8-948f-baceda57a38a
md"""
## Load Graph Data
"""

# ╔═╡ 3676c163-7281-4d20-8923-5f789959db1c
md"Locations of cached Voro-graphs and Crystals with bonding graphs, and the number of graphs to load:"

# ╔═╡ b48da7bd-7daa-4c52-8790-0876558e3e58
md"Array of $nb_graphs structures for which to load graphs:"

# ╔═╡ 3ee4a765-9761-497c-aa20-ddb3099a8ca7
md"""
#### Voro-Graphs
"""

# ╔═╡ 7034dfe5-d22d-4b96-a184-f3c71fd4e515
md"Function to load the Voro-graph for a specific Crystal:"

# ╔═╡ da6a47f9-ea5b-409c-a6e5-c99b7b02f483
md"Array of Voro-graphs for the selected structures:"

# ╔═╡ ddb1e0ed-075d-43fd-a794-b0689a96aedd
md"""
#### Crystals
"""

# ╔═╡ 7fd7397c-6b08-4e1e-92ed-83b63dbb88c2
md"Function to load the Crystal corresponding to a specific Voro-graph:"

# ╔═╡ ab9d16c2-7748-4d17-97aa-61652287efd4
md"Array of Crystals for the selected structures:"

# ╔═╡ 85ddd136-f67b-49c3-b01d-35c594299618
md"""
## Calculate Properties
"""

# ╔═╡ 3d123544-1e2b-488d-a5a3-4caf29684eca
md"""
#### Node Degree
"""

# ╔═╡ da132e30-afe3-4988-9589-0c7bb46a7946
md"""
#### Node Count
"""

# ╔═╡ 23245672-bacc-4c89-b991-7d5fe3620152
node_counts = nv.(voro_graphs);

# ╔═╡ 95baab72-2d16-4534-89ac-972eb0f76d04
md"""
#### Connected Components
"""

# ╔═╡ 43eb7e0d-5ad5-4f93-99de-3234ebde5dc0
conn_comps = length.(connected_components.(voro_graphs));

# ╔═╡ 472d365a-e669-4dd5-8898-1662b70dee30
md"""
#### Graph Diameter
"""

# ╔═╡ 539b09b4-2e43-4f1f-987d-a41afca0d6d6
diameters = filter(d -> d ≠ Inf, [diameter(g) for g in voro_graphs if nv(g) > 0]);

# ╔═╡ dc20dac0-470b-458c-b199-74ad4d69f36a
md"""
#### Accessible Volume Coverage
"""

# ╔═╡ 26fe3e15-50f0-4e77-b4d0-ad273eeccd94
md"The probe molecule and forcefield to use for calculating pore volume:"

# ╔═╡ 5780fdc6-9eb1-4040-aa36-b7f4e2a51d09
begin
	probe = Molecule("He")
	ljforcefield = LJForceField("UFF")
end;

# ╔═╡ 477b609c-3b4a-4262-a6cf-0b5314d8794e
md"A function to calculate the pore volume and accessibility grid of a crystal:"

# ╔═╡ c8e495f8-e181-4ae9-82c3-6279e13f06a8
function pore_grid(crystal::Crystal)
	xtal = deepcopy(crystal)
	remove_bonds!(xtal)
	grid, _, _ = compute_accessibility_grid(xtal, probe, ljforcefield)
	return grid
end;

# ╔═╡ 9ae14fca-e430-416b-992a-9b25428f8abc
md"Arrays of pore volumes and accessibility grids for the selected structures:"

# ╔═╡ d341f92f-00eb-401a-80d7-0f784f57d061
grids = pore_grid.(crystals);

# ╔═╡ f5636586-b30a-4787-abbb-53cb41687049
md"A function to determine if point $(i,j,k)$ of the grid is within the sphere of radius $r$ centered at $x$:"

# ╔═╡ 7b7fa199-c78f-4e89-9375-50364ca3331a
function covered(xtal, accessibility_grid, i, j, k, x, r)
	p = id_to_xf(accessibility_grid.n_pts, (i, j, k))
	d_xp = MLMolGraph.pbc_distance(x, p, xtal.box)
	return d_xp < r
end;

# ╔═╡ a7e908d8-9955-4977-adba-d71a70d05cea
md"A function to calculate the proportion of a Crystal's accessible pore volume covered by the union of Voro-graph spheres:"

# ╔═╡ 4b0807ac-2bcb-47af-bc64-548c499b274c
function sphere_coverage(accessibility_grid, xtal, voro_graph)
	voro_grid = falses(size(accessibility_grid.data))
	n_pts_i, n_pts_j, n_pts_k = accessibility_grid.n_pts
	for v in vertices(voro_graph)
		x = MLMolGraph.shift_back(get_prop(voro_graph, v, :point).coords, xtal.box)
		r = get_prop(voro_graph, v, :radius)
		for k in 1:n_pts_k
			for j in 1:n_pts_j
				for i in 1:n_pts_i
					voro_grid[i, j, k] = 
						voro_grid[i, j, k] || 
						covered(xtal, accessibility_grid, i, j, k, x, r)
				end
			end
		end
	end
	bitvec1 = reshape(accessibility_grid.data .== true, :)
	bitvec2 = reshape(voro_grid, :)
	bitvec3 = bitvec1 .& bitvec2
	return sum(bitvec3) / sum(bitvec1)
end;

# ╔═╡ 1a633440-be96-4396-a4b0-1455857e623f
md"Array of sphere coverages as proportion of each structure's accessible pore volume:"

# ╔═╡ 3b8a77ab-e7b9-4003-8961-412ebcbaafe7
sphere_coverages = sphere_coverage.(grids, crystals, voro_graphs);

# ╔═╡ 644dea49-535e-452f-a4a0-7e22e1b25a83
md"""
## Plot Results
"""

# ╔═╡ 2b64bef1-20cf-451d-99cf-1f3a01cb4cae
begin
	local c = count(d -> d == 0, degree_population)
	md"""
	Graphs w/ nodes of degree zero: $c
	!!! note 
		No Voro-graphs end up w/ nodes having degree zero.  Some graphs have extremely highly-connected nodes (probably not an issue...)
	"""
end

# ╔═╡ 9c834c76-7d0e-4eb2-bb2a-732f0c668903
begin
	#using PyPlot
	figure()
	title("Voro-Graph Node Count Distribution")
	xlabel("Number of Vertices")
	ylabel("Count")
	hist(node_counts, bins=maximum(node_counts)+1)
	gcf()
end

# ╔═╡ a4ae5b16-795d-4d67-9091-eb7a28479351
begin
	local c = count(n -> n == 0, node_counts)
	local p = round(100*c/length(node_counts), digits=2)
md"""
 $(p) % of Voro-graphs have zero nodes.

!!! note
	Some (but not many) Voro-graphs contain zero nodes. Why? Are these graphs w/o He-accessible pore space?  Also, some graphs contain very large numbers of nodes (hopefully not a problem...)
"""
end

# ╔═╡ 6922648d-7cc2-4f65-bb6c-4f3bbf5ed5f2
begin
	figure()
	title("Voro-Graph Connectivity Distribution")
	xlabel("Number of Connected Components")
	ylabel("Count")
	hist(conn_comps, bins=maximum(conn_comps)+1)
	gcf()
end

# ╔═╡ a4bfbbf2-0e5b-4277-aaa8-5126e7266ec9
begin
	local c = count(c -> c == 1, conn_comps)
	local p = 100*round(c/length(conn_comps), digits=2)
md"""
 $p % of graphs are a single connected component.

!!! note
	Very few Voro-graphs contain disjoint subgraphs.  Are these structures with highly distinct 1D channels?
"""
end

# ╔═╡ 5f14a579-beea-4bc1-8908-426eaabf04a1
begin
	figure()
	title("Voro-Graph Eccentricity Distribution")
	xlabel("Graph Diameter")
	ylabel("Count")
	hist(diameters, bins=Int(maximum(diameters)))
	gcf()
end

# ╔═╡ e490a47d-49de-48cc-b245-cafdc4cef4fb
begin
	local m = Int(maximum(diameters))
	local c = count(d -> d > 4, diameters)
	local p = 100*round(c/length(diameters), digits=2)
md"""
Largest graph diameter = $m

 $p % of Voro-graphs have diameter > 4

!!! note
	The diameters of the Voro-graphs are rarely larger than 4. Will 5 message-passing steps be good?
"""
end

# ╔═╡ c262f198-47b3-4c1e-b3ba-1909c4905a3d
begin
	figure()
	title("Pore Volume Coverage Distribution")
	xlabel("Proportion of Coverage")
	ylabel("Count")
	hist(sphere_coverages, bins=20)
	gcf()
end

# ╔═╡ 38ebc164-d519-4d97-b69d-a4b05613b9e7
begin
	local c = count(s -> s < 0.66, sphere_coverages)
	local p = 100 * round(c / length(sphere_coverages), digits=2)
md"""
 $p % of Voro-graphs provide less than 66 % coverage of the accessible pore volume.

!!! note
	We should probably have a 3rd greedy selection stage that adds additional vertices until some threshold is reached.
"""
end

# ╔═╡ Cell order:
# ╠═9ee1c6ad-8660-422c-a1f5-d2b6b93875b3
# ╠═1e19d039-b180-409b-af52-610bec0fc169
# ╟─028398e0-30a4-4ff8-948f-baceda57a38a
# ╟─3676c163-7281-4d20-8923-5f789959db1c
# ╠═ea41da3c-53a4-4d60-8aa1-606903c3c3ae
# ╟─b48da7bd-7daa-4c52-8790-0876558e3e58
# ╠═6f89aa61-cb1c-4c55-b01f-1099e9a59376
# ╟─3ee4a765-9761-497c-aa20-ddb3099a8ca7
# ╟─7034dfe5-d22d-4b96-a184-f3c71fd4e515
# ╠═04970855-74ce-47fc-9642-0edf4563c9ca
# ╟─da6a47f9-ea5b-409c-a6e5-c99b7b02f483
# ╠═0239bc4c-741c-4176-a66c-526398d0fd44
# ╟─ddb1e0ed-075d-43fd-a794-b0689a96aedd
# ╟─7fd7397c-6b08-4e1e-92ed-83b63dbb88c2
# ╠═f1f0f0f9-8e8d-445b-aa53-cb9ba1354ad3
# ╟─ab9d16c2-7748-4d17-97aa-61652287efd4
# ╠═8b1d1c34-46fc-4e75-8a1a-aeb3f19be118
# ╟─85ddd136-f67b-49c3-b01d-35c594299618
# ╟─3d123544-1e2b-488d-a5a3-4caf29684eca
# ╠═cd76fb5e-0a52-4655-8ff2-a1ec580cb68d
# ╟─da132e30-afe3-4988-9589-0c7bb46a7946
# ╠═23245672-bacc-4c89-b991-7d5fe3620152
# ╟─95baab72-2d16-4534-89ac-972eb0f76d04
# ╠═43eb7e0d-5ad5-4f93-99de-3234ebde5dc0
# ╟─472d365a-e669-4dd5-8898-1662b70dee30
# ╠═539b09b4-2e43-4f1f-987d-a41afca0d6d6
# ╟─dc20dac0-470b-458c-b199-74ad4d69f36a
# ╟─26fe3e15-50f0-4e77-b4d0-ad273eeccd94
# ╠═5780fdc6-9eb1-4040-aa36-b7f4e2a51d09
# ╟─477b609c-3b4a-4262-a6cf-0b5314d8794e
# ╠═c8e495f8-e181-4ae9-82c3-6279e13f06a8
# ╟─9ae14fca-e430-416b-992a-9b25428f8abc
# ╠═d341f92f-00eb-401a-80d7-0f784f57d061
# ╟─f5636586-b30a-4787-abbb-53cb41687049
# ╠═7b7fa199-c78f-4e89-9375-50364ca3331a
# ╟─a7e908d8-9955-4977-adba-d71a70d05cea
# ╠═4b0807ac-2bcb-47af-bc64-548c499b274c
# ╟─1a633440-be96-4396-a4b0-1455857e623f
# ╠═3b8a77ab-e7b9-4003-8961-412ebcbaafe7
# ╟─644dea49-535e-452f-a4a0-7e22e1b25a83
# ╟─da66a716-cf5e-4d51-ace2-ea3c671cefa0
# ╟─2b64bef1-20cf-451d-99cf-1f3a01cb4cae
# ╟─9c834c76-7d0e-4eb2-bb2a-732f0c668903
# ╟─a4ae5b16-795d-4d67-9091-eb7a28479351
# ╟─6922648d-7cc2-4f65-bb6c-4f3bbf5ed5f2
# ╟─a4bfbbf2-0e5b-4277-aaa8-5126e7266ec9
# ╟─5f14a579-beea-4bc1-8908-426eaabf04a1
# ╟─e490a47d-49de-48cc-b245-cafdc4cef4fb
# ╟─c262f198-47b3-4c1e-b3ba-1909c4905a3d
# ╟─38ebc164-d519-4d97-b69d-a4b05613b9e7
