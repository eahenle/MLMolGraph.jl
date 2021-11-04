### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 9ee1c6ad-8660-422c-a1f5-d2b6b93875b3
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 1e19d039-b180-409b-af52-610bec0fc169
using Graphs, JLD2, MetaGraphs, MLMolGraph

# ╔═╡ ea41da3c-53a4-4d60-8aa1-606903c3c3ae
dir = "data/cache/vspn"

# ╔═╡ 04970855-74ce-47fc-9642-0edf4563c9ca
function load_graph(jld, dir)
    @load joinpath(dir, jld) obj
    voro_graph, _ = obj
    return voro_graph
end

# ╔═╡ 0239bc4c-741c-4176-a66c-526398d0fd44
graphs = [load_graph(jld, dir) for jld in readdir(dir)]

# ╔═╡ cd76fb5e-0a52-4655-8ff2-a1ec580cb68d
degree_arrays = [degree(g) for g in graphs]

# ╔═╡ 644dea49-535e-452f-a4a0-7e22e1b25a83
degree_population = reduce(vcat, degree_arrays)

# ╔═╡ 048edeb1-5595-4e5f-b1da-012c8b335a24
begin
	using CairoMakie
	fig = Figure()
	ax = Axis(fig[1,1])
	hist!(degree_population)
	fig
end

# ╔═╡ da66a716-cf5e-4d51-ace2-ea3c671cefa0
begin
	using PyPlot
	figure()
	title("Degree Distribution")
	xlabel("Node Degree")
	ylabel("Frequency")
	hist(degree_population, density=true, bins=250)
	gcf()
end

# ╔═╡ 2b64bef1-20cf-451d-99cf-1f3a01cb4cae
count(d -> d == 0, degree_population)

# ╔═╡ 1df8a8a6-813a-4cdd-9d6f-c8fbba83292f
md"""
!!! note 
	No Voro-graphs end up w/ nodes having degree zero.  Some graphs have extremely highly-connected nodes (probably not an issue...)
"""

# ╔═╡ 8d05826d-e6a4-4192-af3d-068b84fbf053
node_counts = [nv(g) for g in graphs]

# ╔═╡ 25950ceb-18ab-4f65-88dd-7e8dcc016dc6
sum(node_counts)

# ╔═╡ 9c834c76-7d0e-4eb2-bb2a-732f0c668903
begin
	#using PyPlot
	figure()
	title("Vertex Count Distribution")
	xlabel("Number of Vertices")
	ylabel("Frequency")
	hist(node_counts, density=true, bins=600)
	gcf()
end

# ╔═╡ 2d58d822-6163-4454-a71f-58b6c73923a2
count(n -> n == 0, node_counts)

# ╔═╡ a4ae5b16-795d-4d67-9091-eb7a28479351
md"""
!!! note
	Some Voro-graphs contain zero nodes! Why? Are these graphs w/o He-accessible pore space?  Also, some graphs contain very large numbers of nodes (hopefully not a problem...)
"""

# ╔═╡ Cell order:
# ╠═9ee1c6ad-8660-422c-a1f5-d2b6b93875b3
# ╠═1e19d039-b180-409b-af52-610bec0fc169
# ╠═ea41da3c-53a4-4d60-8aa1-606903c3c3ae
# ╠═04970855-74ce-47fc-9642-0edf4563c9ca
# ╠═0239bc4c-741c-4176-a66c-526398d0fd44
# ╠═cd76fb5e-0a52-4655-8ff2-a1ec580cb68d
# ╠═644dea49-535e-452f-a4a0-7e22e1b25a83
# ╠═048edeb1-5595-4e5f-b1da-012c8b335a24
# ╠═da66a716-cf5e-4d51-ace2-ea3c671cefa0
# ╠═2b64bef1-20cf-451d-99cf-1f3a01cb4cae
# ╟─1df8a8a6-813a-4cdd-9d6f-c8fbba83292f
# ╠═8d05826d-e6a4-4192-af3d-068b84fbf053
# ╠═25950ceb-18ab-4f65-88dd-7e8dcc016dc6
# ╠═9c834c76-7d0e-4eb2-bb2a-732f0c668903
# ╠═2d58d822-6163-4454-a71f-58b6c73923a2
# ╟─a4ae5b16-795d-4d67-9091-eb7a28479351
