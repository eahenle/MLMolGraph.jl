function xtals2primitive(xtal_list::Vector{String})
    xtal_paths = [joinpath(rc[:cache][:primitive], xtal) for xtal ∈ xtal_list]
    @sync @distributed for xtal_file ∈ xtal_paths
        cached(xtal_file) do 
            try
                xtal = Crystal(xtal_file, remove_duplicates=true)
                return primitive_cell(xtal)
            catch
                return Crystal("bad input", unit_cube(), 
                    Atoms([:foo], Frac([0.;0.;0.])), 
                    Charges{Frac}(0))
            end
        end
    end    
end


function isgood(xtal_file::String)::Bool
    @load joinpath(rc[:cache][:bonded_xtals], xtal_file) obj
    _, good = obj
    return good
end


function _bondNclassify(xtal_name::String, primitive_cache::String, bonded_cache::String)::Tuple{Crystal,Bool}
    cached(joinpath(bonded_cache, xtal_name)) do 
        @load joinpath(primitive_cache, xtal_name) obj
        xtal = obj
        if infer_bonds!(xtal, true)
            good = true
        else
            good = false
        end
        return (xtal, good)
    end
end


function bondNclassify(xtal_list::Vector{String})::Vector{String}
    l = length(xtal_list)
    pcs = [rc[:cache][:primitive] for _ ∈ 1:l]
    bcs = [rc[:cache][:bonded_xtals] for _ ∈ 1:l]
    pmap(_bondNclassify, xtal_list, pcs, bcs)
    return [xtal_file for xtal_file ∈ xtal_list if isgood(xtal_file)]
end


function encode(xtal_list::Vector{String})::Tuple{Dict{Symbol,Int},Int}
	element_to_int = Dict{Symbol,Int}()
	nb_elements = 0
    max_valency = 0
	for xtal_file ∈ xtal_list
		@load joinpath(rc[:cache][:bonded_xtals], xtal_file) obj
        xtal, good = obj
        if good
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
	end
	return element_to_int, max_valency
end


function read_targets(source::String, examples::Vector{String}, target_symbol::Symbol)::DataFrame
    df = CSV.read(joinpath(rc[:paths][:data], source), DataFrame, delim=", ")
    select!(df, [:name, target_symbol])
    filter!(r -> "$(r.name).cif" ∈ examples, df)
    return df
end


function node_feature_matrix(xtal::Crystal, max_valency::Int, element_to_int::Dict{Symbol,Int})
    embedding_length = length(element_to_int) + max_valency
    X = zeros(Int, xtal.atoms.n, embedding_length)
    for (i, atom) in enumerate(xtal.atoms.species)
        X[i, element_to_int[atom]] = 1
        X[i, embedding_length - degree(xtal.bonds)[i] + 1] = 1
    end
    return X
end


function edge_vectors(graph::MetaGraph, args::Dict{Symbol,Any})::Union{Tuple{Vector{Int},Vector{Int},Vector{Float64}},Tuple{Vector{Int},Vector{Int},Vector{Int},Vector{Float64}}}
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
    if args[:bonds]
        return edg_srcs, edg_dsts, edg_lens
    elseif args[:vspn]
        edg_type = zeros(Int, 1, l)
        edg_wght = zeros(Float64, 1, l)
        for (i, edge) in enumerate(edges(graph))
            s = src(edge)
            d = dst(edge)
            type = get_prop(graph, s, d, :type)
            if type == :AA
                edg_type[i] = edg_type[i + edge_count] = 1
                edg_wght[i] = edg_wght[i + edge_count] = get_prop(graph, s, d, :distance)
            elseif type == :AV
                edg_type[i] = edg_type[i + edge_count] = 2
                edg_wght[i] = edg_wght[i + edge_count] = get_prop(graph, s, d, :distance)
            else # type == :VV
                edg_type[i] = edg_type[i + edge_count] = 3
                edg_wght[i] = edg_wght[i + edge_count] = get_prop(graph, s, d, :lens)
            end
        end
        return edg_srcs, edg_dsts, edg_type, edg_wght
    end
end


function bond_angle_vecs(xtal::Crystal)::Tuple{Vector{Int},Vector{Int},Vector{Int},Vector{Float64}}
    I = Int[]
    J = Int[]
    K = Int[]
    θ = Float64[]
    for BA ∈ edges(xtal.bonds)
        A = dst(BA)
        B = src(BA)
        for C ∈ neighbors(xtal.bonds, B)
            if A == C # only need actual 3-body angles
                continue
            end
            α = bond_angle(xtal, A, B, C)
            # ∠ABC = ∠CBA
            append!(I, A)
            append!(J, B)
            append!(K, C)
            append!(θ, α)
            append!(I, C)
            append!(J, B)
            append!(K, A)
            append!(θ, α)
        end
    end
    return I, J, K, θ
end


function write_data(xtal::Crystal, name::String, element_to_int::Dict{Symbol,Int}, max_valency::Int, graphs_path::String, args::Dict{Symbol,Any}, config::Union{Nothing,VSPNConfig}=nothing)
    X = node_feature_matrix(xtal, max_valency, element_to_int)
    X_name = chop(name, tail=4)
	# bond graph
    if args[:bonds]
        A, B, D = edge_vectors(xtal.bonds, args)
        npzwrite(joinpath(graphs_path, X_name * "_node_features.npy"), X)
        npzwrite(joinpath(graphs_path, X_name * "_edges_src.npy"), A)
        npzwrite(joinpath(graphs_path, X_name * "_edges_dst.npy"), B)
        npzwrite(joinpath(graphs_path, X_name * "_euc.npy"), D)
    end
    # VSPN graph
    if args[:vspn]
        vspn = vspn_graph(xtal, config)
        A, B, T, W = edge_vectors(vspn, args)
        npzwrite(joinpath(rc[:paths][:graphs], X_name * "_vspn_node_features.npy"), X)
        npzwrite(joinpath(rc[:paths][:graphs], X_name * "_vspn_edges_src.npy"), A)
        npzwrite(joinpath(rc[:paths][:graphs], X_name * "_vspn_edges_dst.npy"), B)
        npzwrite(joinpath(rc[:paths][:graphs], X_name * "_vspn_edges_types.npy"), T)
        npzwrite(joinpath(rc[:paths][:graphs], X_name * "_vspn_edges_weights.npy"), W)
    end
	# bond angles
    if args[:angles]
        I, J, K, θ = bond_angle_vecs(xtal)
        npzwrite(joinpath(graphs_path, X_name * "_angles_I.npy"), I)
        npzwrite(joinpath(graphs_path, X_name * "_angles_J.npy"), J)
        npzwrite(joinpath(graphs_path, X_name * "_angles_K.npy"), K)
        npzwrite(joinpath(graphs_path, X_name * "_angles_theta.npy"), θ)
    end
end


function process_example(xtal_name::String, element_to_int::Dict{Symbol,Int}, max_valency::Int, bonded_xtals_cache::String, graphs_path::String, args::Dict{Symbol,Any}, config::Union{Nothing,VSPNConfig}=nothing)
    @load joinpath(bonded_xtals_cache, xtal_name) obj
    xtal, _ = obj
    write_data(xtal, xtal_name, element_to_int, max_valency, graphs_path, args, config=nothing)
end


function process_examples(good_xtals::Vector{String}, element_to_int::Dict{Symbol,Int}, max_valency::Int, args::Dict{Symbol,Any})
    l = length(good_xtals)
    els = [element_to_int for _ ∈ 1:l]
    mvs = [max_valency for _ ∈ 1:l]
    bxc = [rc[:cache][:bonded_xtals] for _ ∈ 1:l]
    gps = [rc[:paths][:graphs] for _ ∈ 1:l]
    ags = [args for _ ∈ 1:l]
    if args[:vspn]
        config = VSPNConfig(args[:probe], args[:forcefield])
        cfs = [config for _ ∈ 1:l]
        pmap(process_example, good_xtals, els, mvs, bxc, gps, ags, cfs)
    else
        pmap(process_example, good_xtals, els, mvs, bxc, gps, ags)
    end
end