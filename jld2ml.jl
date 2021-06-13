using Distributed

@everywhere using CSV, DataFrames, JLD2, LightGraphs, MetaGraphs, NPZ, Xtals

# settings and paths
@everywhere begin
    data_root = joinpath(pwd(), "Mercado_COFs") # where to find data
    target_col_name = "deliverable capacity [v STP/v]" # target for model
    target_inputfile = "properties.csv" # name of file with target data
    bondcache = joinpath(data_root, "bonded_xtal_jlds")  # location of JLD2 cache for bonded crystals
end


# return a dataframe w/ crystal names and target values for all materials with good bonds
function read_targets(target_symbol::Symbol)::DataFrame
    good_xtals = readdir(bondcache)
    df = CSV.read(joinpath(data_root, target_inputfile), DataFrame, delim=", ")#[1:10, :]
    select!(df, [:name, target_symbol])
    filter!(r -> "$(r.name).cif.jld2" ∈ good_xtals, df)
    return df
end


# create dictionary for the one-hot encoding column of atom species
function build_encoding_data(df::DataFrame)::Tuple{Dict{Symbol,Int},Int}
	element_to_int = Dict{Symbol,Int}()
	nb_elements = 0
    max_valency = 0
	for material ∈ df.name
		@load joinpath(bondcache, "$material.cif.jld2") xtal
		unique_atoms = unique(xtal.atoms.species)
		for (i,atom) ∈ enumerate(unique_atoms)
            n = degree(xtal.bonds)[i]
            if n > max_valency
                max_valency = n
            end
			if atom ∉ keys(element_to_int)
				nb_elements += 1
				element_to_int[atom] = nb_elements
			end
		end
	end
	return element_to_int, max_valency
end


# node feature matrix. cols: 1-hot encoding of element type and 1-hot of valency (1-6)
# nodes given by [num_nodes, num_features] array (each node has num_features features)
@everywhere function node_feature_matrix(xtal::Crystal, max_valency::Int)::Matrix{Int}
    embedding_length = length(element_to_int) + max_valency
    X = zeros(Int, xtal.atoms.n, embedding_length)
    for (i, atom) in enumerate(xtal.atoms.species)
        X[i, element_to_int[atom]] = 1
        X[i, embedding_length - degree(xtal.bonds)[i] + 1] = 1
    end
    return X
end


# first vector: edge source nodes
# second vector: edge destination nodes
# third vector: distance between bonded nodes
@everywhere function edge_vectors(graph::AbstractGraph)
    edge_count = ne(graph)
	l = 2 * edge_count
    edg_srcs = zeros(Int, 1, l)
    edg_dsts = zeros(Int, 1, l)
    edg_lens = zeros(Float64, 1, l)
    for (i, edge) in enumerate(edges(graph))
        edg_srcs[i] = edg_dsts[i + edge_count] = edge.src - 1 # python numbering
        edg_dsts[i] = edg_srcs[i + edge_count] = edge.dst - 1
        edg_lens[i] = edg_lens[i + edge_count] = get_prop(graph, edge, :distance)
    end
    return edg_srcs, edg_dsts, edg_lens
end


@everywhere function bond_angle_matrix(xtal::Crystal)
    angles = Array{Float64}(undef, xtal.atoms.n, xtal.atoms.n, xtal.atoms.n)
    for BA ∈ edges(xtal.bonds)
        A = dst(BA)
        B = src(BA)
        for BC ∈ edges(xtal.bonds)
            C = dst(BC)
            angles[A, B, C] = angles[C, B, A] = bond_angle(xtal, A, B, C)
        end
    end
    return angles
end


# record ML inputs in numpy format
@everywhere function write_data(xtal::Crystal, max_valency::Int)
    X = node_feature_matrix(xtal, max_valency)
    X_name = chop(xtal.name, tail=4)
	# bond graph
	A, B, D = edge_vectors(xtal.bonds)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_node_features.npy"), X)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_edges_src.npy"), A)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_edges_dst.npy"), B)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_euc.npy"), D)
	# bond angles
    angleABC = bond_angle_matrix(xtal)
    npzwrite(joinpath(pwd(), "Graphs", X_name * "_angles.npy"), angleABC)
end


# main block
@everywhere begin
    target_symbol = Symbol(target_col_name)
    df = read_targets(target_symbol)
    element_to_int, max_valency = build_encoding_data(df)
end
@info "Input data" df
@info "Encoding" element_to_int max_valency
@distributed for material in readdir(bondcache)
    @load joinpath(bondcache, material) xtal
    write_data(xtal, max_valency)
end
npzwrite("y.npy", df[:, target_symbol])
CSV.write(joinpath("atom_to_int.csv"), element_to_int)
CSV.write(joinpath("cofs.csv"), df)
