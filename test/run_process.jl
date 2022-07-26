module Run_process_test

using CSV, DataFrames, Graphs, MetaGraphs, PyCall, Test, MLMolGraph, Xtals, XtalsPyTools

target = "working_capacity_vacuum_swing [mmol/g]"

temp_dir = joinpath(pwd(), ".temp")
if !isdir(temp_dir)
    mkpath(temp_dir)
end

"""
    g = rebuild_graph(xn)

Reconstruct a graph (referenced by structure name) from ML input data
"""
function rebuild_graph(xtal_name::String, int_to_atom::Dict{Int, Symbol}, temp_dir::String)::MetaGraph
    # read ML data files
	xtal_data = MLMolGraph.load_data(xtal_name, temp_dir)
    X = xtal_data[:atom_features]
	es = xtal_data[:bond_edges][1]
	ed = xtal_data[:bond_edges][2]
    @assert length(es) == length(ed) && length(es) > 0 "length mismatch in edge arrays for $xtal_name"
    # reconstruct the graph
	g = MetaGraph(size(X, 1))
    # copy node labels
    for (i, row) in enumerate(eachrow(X))
        @assert set_prop!(g, i, :species, int_to_atom[findfirst(row .== 1)])
    end
    # copy edges
	for (s, d) in zip(es, ed)
        if !has_edge(g, d, s)
            @assert s <= size(X, 1) && d <= size(X, 1) && s >= 1 && d >= 1 "invalid index in edge arrays for $xtal_name: s = $s ($(maximum(es))), d = $d ($(maximum(ed))), [X] = $(size(X))"
            @assert add_edge!(g, s, d) nv(g), s, d
        end
	end
	return g
end


@testset "bonds" begin
    # set up args
    arg_dict = make_arg_dict(target)
    arg_dict[:bonds] = true

    # process the crystals
    good_xtals = run_process(arg_dict, temp_dir=temp_dir)
    # test that all inputs return usable results
    @test length(good_xtals) == 50

    # load crystals
    xtals = Crystal.(good_xtals .* [".cif"])
    # load feature matrices
    py"""
	import pickle
	file = open("_DATASET_.PKL", "rb")
	dataset = pickle.load(file)
	file.close()
	"""
	dataset = py"dataset"
    atom_x = dataset["atom_x"]
    Xs = [[[get(get(atom_x[graph], node-1), i-1).item() for i in 1:length(get(atom_x[graph], node-1))] for node in 1:length(atom_x[graph])] for graph in eachindex(atom_x)]

    # test that feature matrix rows are one-hot encodings
    @test all([all([sum(Xs[i][j]) .== 1 for j in eachindex(Xs[i])]) for i in eachindex(Xs)])

    # convert to primitive cells (to keep consistent w/ processing)
    xtals = primitive_cell.(xtals)

    # read atom encoding dictionary
    atom_to_int = Dict([Symbol(value) => key for (value, key) in dataset["atom_to_int"]])
    # test that atom encodings match original species by atom_to_int mapping for all structures
    @test all([all([atom_to_int[xtals[i].atoms.species[j]] == findfirst(Xs[i][j] .== 1) for j in eachindex(Xs[i])]) for i in eachindex(Xs)])

    # build bonding graphs
    good_bonds = infer_bonds!.(xtals, [true])
    @assert all(good_bonds) # make sure the bonding went right
    # invert encoding dictionary
    int_to_atom = Dict([label => key for (key, label) in atom_to_int])
    @assert length(int_to_atom) == length(atom_to_int) # check bijectivity

    # reconstruct graphs from numpy data
    graphs = rebuild_graph.(good_xtals, [int_to_atom], temp_dir)
    @test all([xtal.bonds for xtal in xtals] .== graphs) # test reconstructed graphs against originals

    # check target data against source CSV
    properties_df = CSV.read(joinpath(rc[:paths][:data], "properties.csv"), DataFrame)
    targets = [y.item() for y in dataset["y"]]
    test_flags = falses(length(properties_df[:, target]))
    for (i, name) in enumerate(good_xtals)
        # extract values from independent locations by name
        pt = filter(row -> row.name == name, properties_df)[:, target][1]
        tt = targets[findfirst(dataset["xtal_name"] .== name)]
        # confirm that targets in each location are the same
        test_flags[i] = isapprox(pt, tt, atol=1e-4)
    end
    @test all(test_flags)
end

@testset "element fractions" begin
    py"""
	import pickle
	file = open("_DATASET_.PKL", "rb")
	dataset = pickle.load(file)
	file.close()
	"""
	dataset = py"dataset"
    @test all(sum.(fs for fs in dataset["elemental_fractions"]) .≈ 1)
end

end
