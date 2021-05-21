### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ d241db54-2c70-11eb-09a8-c7941c678a82
begin
	using Xtals, LightGraphs, DataFrames, CSV, NPZ, JLD2, MetaGraphs
end

# ╔═╡ e9d041e8-2c70-11eb-23e1-858ad08986e4
md"
source of data: [Materials Cloud](https://archive.materialscloud.org/record/2018.0003/v2) [v2]
* `structures.tgz`: crystal structures
* `properties.tgz`: properties

from reference: [Mercado et al.](https://pubs.acs.org/doi/pdf/10.1021/acs.chemmater.8b01425)
"

# ╔═╡ 3d84b88c-2c71-11eb-0d38-5f94ac70c983
begin
	# tell Xtals.jl where to find the .cif files
	set_path_to_crystals(joinpath(pwd(), "structures"))
	if !isdir("Graphs")
		mkdir("Graphs")
	end
	# for caching
	# needs initial null string (for jld2 base dir, aka pwd)
	jldirs = ["", "bonds", "skeleton", "xtal", "bad_bonds"]
	for dir in jldirs
		dir = "jld2/" * dir
		if !isdir(dir)
			mkdir(dir)
		end
	end
end

# ╔═╡ 911ca9f8-2c71-11eb-0601-d5de2c50666d
md"read in properties"

# ╔═╡ 662539ba-2c71-11eb-13de-191a8446da4b
df = CSV.read("properties.csv", DataFrame)#[1:100, :] # subset for quick tests

# ╔═╡ 0a5c3b30-5a8c-11eb-1f93-b39d1e4acc59
md"every column name and string value in the dataframe has a leading space!"

# ╔═╡ 0dd35dbe-2c7b-11eb-27c2-75e563191db9
target_col_name = Symbol(" deliverable capacity [v STP/v]")

# ╔═╡ 83bf5860-582d-11eb-0fe9-f14f9b2892dd
select!(df, [Symbol(" name"), target_col_name])

# ╔═╡ 72ad0ac0-5a87-11eb-27d9-3bf941a4358b
for i ∈ 1:nrow(df)
	df[i, Symbol(" name")] = split(df[i, Symbol(" name")], " ")[2]
end

# ╔═╡ 36535630-5a8a-11eb-132a-79e1981784e4
df

# ╔═╡ b79b4bc0-582e-11eb-3dc0-8133b1f44b54
md"get crystal data, with auto-correcting"

# ╔═╡ d83b52f0-5624-11eb-0b67-a3b43a07363b
begin
	bad_struct = []
	for (m, material) in enumerate(df[:, Symbol(" name")])
		if m % round(nrow(df)/100) == 0 # "progress meter"
			println("$m / $(nrow(df))")
		end
		try
			# cache auto-corrected Crystals
			if !isfile("jld2/xtal/$material.jld2")
				xtal = Crystal(material * ".cif", remove_duplicates=true)
				@save "jld2/xtal/$material.jld2" xtal
			end
		catch exception
			@warn "mat'l #$m ($material) not loaded due to exception" exception
			push!(bad_struct, material)
		end
	end
end

# ╔═╡ 9e252100-5830-11eb-2634-ab3d9026c29f
md"drop materials with bad structures"

# ╔═╡ 36fb5400-582f-11eb-2b50-83f440c80cbf
filter!(row -> !(row[" name"] ∈ bad_struct), df)

# ╔═╡ 9bead0e4-2f66-11eb-230a-351d5969b2c0
md"find unique elements"

# ╔═╡ ec3c3570-5832-11eb-234c-3fae45a58582
function build_element_to_int_mapping()
	element_to_int = Dict{Symbol, Int}()
	nb_elements = 0
	for (m, material) in enumerate(df[:, Symbol(" name")])
		if m % round(nrow(df)/100) == 0 # "progress meter"
			println("$m / $(nrow(df))")
		end
		if !isfile("jld2/xtal/$material.jld2")
			continue
		else
			xtal = nothing
			@load "jld2/xtal/$material.jld2" xtal
		end
		# id all unique atoms in the data set
		unique_atoms = unique(xtal.atoms.species)
		for ua in unique_atoms
			if ! (ua in keys(element_to_int))
				nb_elements += 1
				element_to_int[ua] = nb_elements
			end
		end
	end
	return element_to_int
end

# ╔═╡ f2dcfc70-5832-11eb-2042-5d6a37490ab3
if !isfile("jld2/element_to_int.jld2")
	element_to_int =  build_element_to_int_mapping()
	@save "jld2/element_to_int.jld2" element_to_int
else
	element_to_int = nothing
	@load "jld2/element_to_int.jld2" element_to_int
end
# element_to_int[:C] gives integer for C atoms

# ╔═╡ 16a4b9d0-5aa5-11eb-2467-b1166d5ec119
element_to_int

# ╔═╡ e8fe1b9c-2f53-11eb-3370-b109901239d7
md"
$f(A, X) \rightarrow y\in \mathbb{R}$
"

# ╔═╡ 20e442b4-3c1c-11eb-19ef-eddecda42ab5
# node feature matrix. cols: 1-hot encoding of element type and 1-hot of valency (1-6)
# nodes given by [num_nodes, num_features] array (each node has num_features features)
function node_feature_matrix(xtal::Crystal)
    embedding_length = length(element_to_int) + 4 # for valency
    X = zeros(Int, xtal.atoms.n, embedding_length)
    for (i, atom) in enumerate(xtal.atoms.species)
        X[i, element_to_int[atom]] = 1
        X[i, embedding_length - degree(xtal.bonds)[i] + 1] = 1
    end
    return X
end

# ╔═╡ efe83000-5aaf-11eb-2a91-37c13676c15a
function edge_vectors(xtal::Crystal, graph::AbstractGraph)
    edge_count = ne(graph)
	l = 2 * edge_count
    edg_srcs = zeros(Int, 1, l)
    edg_dsts = zeros(Int, 1, l)
    edg_lens = zeros(Float64, 1, l)
    for (i, edge) in enumerate(edges(graph))
        edg_srcs[i] = edg_dsts[i + edge_count] = edge.src - 1
        edg_dsts[i] = edg_srcs[i + edge_count] = edge.dst - 1
        edg_lens[i] = edg_lens[i + edge_count] = get_prop(graph, edge, :distance)
			#distance(xtal.atoms, xtal.box, edge.src, edge.dst, true)
    end
    return edg_srcs, edg_dsts, edg_lens
end

# ╔═╡ 7a14cea0-5ab0-11eb-1a98-7f49ba6c1fc5
edge_vectors(xtal::Crystal) = edge_vectors(xtal, xtal.bonds)

# ╔═╡ 2ee27fb6-44d7-11eb-0df7-adfb8f1584d4
md"build bonds"

# ╔═╡ 77ad0cc8-2c79-11eb-12c1-bf8d9161d476
begin
	# cached bond building
	for (m, material) in enumerate(eachrow(df))
		if m % round(nrow(df)/100) == 0 # "progress meter"
			println("$m / $(nrow(df))")
		end
		if !isfile("jld2/bonds/$(material[Symbol(" name")]).jld2")
			if !isfile("jld2/xtal/$(material[Symbol(" name")]).jld2")
				continue
			else
				xtal = nothing
				@load "jld2/xtal/$(material[Symbol(" name")]).jld2" xtal
			end
			infer_bonds!(xtal, true)
			@save "jld2/bonds/$(material[Symbol(" name")]).jld2" xtal
		end
	end
end

# ╔═╡ 67eace6e-5833-11eb-1a7e-2b6963857965
md"find bad bonds"

# ╔═╡ 48244c60-5621-11eb-2b53-5f9cf3cae982
begin
	bad_bonds = []
	for matl ∈ df[:, Symbol(" name")]
		xtal = nothing
		@load "jld2/bonds/$matl.jld2" xtal
		if !bond_sanity_check(xtal)
			push!(bad_bonds, matl)
		end
	end
end

# ╔═╡ 6f8d9c20-5833-11eb-2fba-61e1d027a11a
md"drop bad bonds from df"

# ╔═╡ 76a7d8de-5833-11eb-0b67-09807bbd755a
filter!(row -> !(row[:, Symbol(" name")] ∈ bad_bonds), df)

# ╔═╡ 86e038b0-5833-11eb-23f7-7737f6f64003
md"""build the "skeleton" graphs"""

# ╔═╡ 981ce920-5833-11eb-1f77-8706f8487482
function get_skeleton(xtal::Crystal)::MetaGraph
	# copy the nodes
	skeleton = MetaGraph(nv(xtal.bonds))
	V = collect(vertices(xtal.bonds))

	# make edges that are not in bonding graph
	for v ∈ V
		for w ∈ V
			if v == w || has_edge(xtal.bonds, v, w) || has_edge(xtal.bonds, w, v)
				continue
			else
				add_edge!(skeleton, v, w)
				set_prop(xtal.bonds, (v, w), :distance, distance(xtal.atoms, xtal.box, edge.src, edge.dst, true))
			end
		end
	end

	return skeleton
end

# ╔═╡ f2643ce0-5567-11eb-2157-4b3727b15b1e
for (m, matl) in enumerate(df[:, Symbol(" name")])
	if m % round(nrow(df)/100) == 0 # "progress meter"
		println("$m / $(nrow(df))")
	end
	if isfile("jld2/skeleton/$matl.jld2") || !isfile("jld2/bonds/$matl.jld2")
		continue
	end
	xtal = nothing
	@load "jld2/bonds/$matl.jld2" xtal
	skeleton = get_skeleton(xtal)
	@save "jld2/skeleton/$matl.jld2" skeleton
end

# ╔═╡ 3955cc80-5834-11eb-0185-b5e967253fd5
md"prepare ML data"

# ╔═╡ 736c4850-56f3-11eb-2576-dd979e8072ef
function write_data(xtal::Crystal, skeleton::AbstractGraph)
    X = node_feature_matrix(xtal)
    X_name = replace(xtal.name, ".cif" => "")
	# bond graph
	A, B, D = edge_vectors(xtal)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_node_features.npy"), X)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_edges_src.npy"), A)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_edges_dst.npy"), B)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_euc.npy"), D)
	# skeleton graph
	A, B, D = edge_vectors(xtal, skeleton)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_skeleton_edges_src.npy"), A)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_skeleton_edges_dst.npy"), B)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_skeleton_euc.npy"), D)
end

# ╔═╡ 171e06e0-5623-11eb-13b0-319b44bb9a2d
begin
	# target vector
	y = Float64[]
	for (m, material) in enumerate(eachrow(df))
		if m % round(nrow(df)/100) == 0 # "progress meter"
			println("$m / $(nrow(df))")
		end
		xtal = nothing
		skeleton = nothing
		@load "jld2/bonds/" * material[:, Symbol(" name")] * ".jld2" xtal
		@load "jld2/skeleton/" * material[:, Symbol(" name")] * ".jld2" skeleton
		write_data(xtal, skeleton)
		npzwrite()
		push!(y, material[target_col_name])
	end
	npzwrite("y.npy", y) # save target vector
	CSV.write(joinpath("atom_to_int.csv"), element_to_int)
	CSV.write(joinpath("cofs.csv"), df)
end

# ╔═╡ Cell order:
# ╠═d241db54-2c70-11eb-09a8-c7941c678a82
# ╟─e9d041e8-2c70-11eb-23e1-858ad08986e4
# ╠═3d84b88c-2c71-11eb-0d38-5f94ac70c983
# ╟─911ca9f8-2c71-11eb-0601-d5de2c50666d
# ╠═662539ba-2c71-11eb-13de-191a8446da4b
# ╟─0a5c3b30-5a8c-11eb-1f93-b39d1e4acc59
# ╠═0dd35dbe-2c7b-11eb-27c2-75e563191db9
# ╠═83bf5860-582d-11eb-0fe9-f14f9b2892dd
# ╠═72ad0ac0-5a87-11eb-27d9-3bf941a4358b
# ╠═36535630-5a8a-11eb-132a-79e1981784e4
# ╟─b79b4bc0-582e-11eb-3dc0-8133b1f44b54
# ╠═d83b52f0-5624-11eb-0b67-a3b43a07363b
# ╟─9e252100-5830-11eb-2634-ab3d9026c29f
# ╠═36fb5400-582f-11eb-2b50-83f440c80cbf
# ╟─9bead0e4-2f66-11eb-230a-351d5969b2c0
# ╠═ec3c3570-5832-11eb-234c-3fae45a58582
# ╠═f2dcfc70-5832-11eb-2042-5d6a37490ab3
# ╠═16a4b9d0-5aa5-11eb-2467-b1166d5ec119
# ╟─e8fe1b9c-2f53-11eb-3370-b109901239d7
# ╠═20e442b4-3c1c-11eb-19ef-eddecda42ab5
# ╠═efe83000-5aaf-11eb-2a91-37c13676c15a
# ╠═7a14cea0-5ab0-11eb-1a98-7f49ba6c1fc5
# ╟─2ee27fb6-44d7-11eb-0df7-adfb8f1584d4
# ╠═77ad0cc8-2c79-11eb-12c1-bf8d9161d476
# ╟─67eace6e-5833-11eb-1a7e-2b6963857965
# ╠═48244c60-5621-11eb-2b53-5f9cf3cae982
# ╟─6f8d9c20-5833-11eb-2fba-61e1d027a11a
# ╠═76a7d8de-5833-11eb-0b67-09807bbd755a
# ╟─86e038b0-5833-11eb-23f7-7737f6f64003
# ╠═981ce920-5833-11eb-1f77-8706f8487482
# ╠═f2643ce0-5567-11eb-2157-4b3727b15b1e
# ╟─3955cc80-5834-11eb-0185-b5e967253fd5
# ╠═736c4850-56f3-11eb-2576-dd979e8072ef
# ╠═171e06e0-5623-11eb-13b0-319b44bb9a2d
