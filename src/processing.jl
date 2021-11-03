function loadcheck(xtal_list::Vector{String}, tolerance::Float64)
    @showprogress "Checking file loadability:" @distributed for xtal_file ∈ xtal_list
        xtal_filename = String(split(xtal_file, "/")[end])
        cached("primitive/$xtal_filename") do 
            try
                return Crystal(xtal_filename, remove_duplicates=true, net_charge_tol=tolerance)
            catch exception
                @error xtal_file exception
                return Crystal("bad input", unit_cube(), 
                    Atoms([:foo], Frac([0.;0.;0.])), 
                    Charges{Frac}(0))
            end
        end
    end
end


# TODO add optional arguments for changing the net_charge_tol
function xtals2primitive(xtal_list::Vector{String}, tolerance::Float64)
    @showprogress "Converting to primitive cells:" @distributed for xtal_file ∈ xtal_list
        xtal_filename = String(split(xtal_file, "/")[end])
        cached("primitive/$xtal_filename") do 
            try
                xtal = Crystal(xtal_filename, remove_duplicates=true, net_charge_tol=tolerance)
                return primitive_cell(xtal)
            catch exception
                @error xtal_file exception
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


function _bondNclassify(xtal_name::String)::Tuple{Crystal,Bool}
    cached("bonded_xtals/$xtal_name") do 
        @load joinpath(rc[:cache][:primitive], xtal_name) obj
        xtal = obj
        if infer_bonds!(xtal, true, calculate_vectors=true)
            good = true
        else
            good = false
        end
        return (xtal, good)
    end
end


function bondNclassify(xtal_list::Vector{String})::Vector{String}
    @showprogress "Inferring bonds:" pmap(xtal_list) do x
        _bondNclassify(x)
    end
    good = SharedArray{Bool}(length(xtal_list))
    @showprogress "Classifying structures:" @distributed for i ∈ 1:length(xtal_list)
        good[i] = isgood(xtal_list[i])
    end
    return xtal_list[good]
end


function unique_elements(xtal_file::String, args)
	@load joinpath(rc[:cache][:bonded_xtals], xtal_file) obj
    xtal = obj[1]
    elements = unique(xtal.atoms.species)
    if args[:verbose]
        @info "$elements unique elements in $(xtal.name)"
    end
    return elements
end


function encode(xtal_list::Vector{String}, args)
    n = length(xtal_list)
    all_elements = [keys(rc[:covalent_radii])...] # list of all elements known to Xtals
    elements = SharedArray{Bool}(length(all_elements)) # elements[i] == true iff keys(all_elements)[i] ∈ some xtal
    if args[:verbose]
        @info "Encoding $n structures' atoms against list of $(length(all_elements)) elements"
    end
    @showprogress "Determining encoding scheme:" @distributed for i ∈ 1:n
        xtal_elements = unique_elements(xtal_list[i], args)
        elements[[findfirst(isequal(element), all_elements) for element in xtal_elements]] .= true
    end
	element_to_int = Dict{Symbol,Int}([element => i for (i, element) ∈ enumerate(all_elements[elements])])
    encoding_length = length(element_to_int)
	return element_to_int, encoding_length
end


function read_targets(source::String, examples::Vector{String}, target_symbol::Symbol)::DataFrame
    df = CSV.read(joinpath(rc[:paths][:data], source), DataFrame, delim=",")
    select!(df, [:name, target_symbol])
    filter!(r -> "$(r.name).cif" ∈ examples, df)
    return df
end


function node_feature_matrix(xtal::Crystal, element_to_int::Dict{Symbol,Int})
    embedding_length = length(element_to_int)
    X = zeros(Int, xtal.atoms.n, embedding_length)
    for (i, atom) in enumerate(xtal.atoms.species)
        X[i, element_to_int[atom]] = 1
        X[i, embedding_length - degree(xtal.bonds)[i] + 1] = 1
    end
    return X
end


function edge_vectors(graph::MetaGraph, args::Dict{Symbol,Any})
    edge_count = ne(graph)
	l = 2 * edge_count
    edg_srcs = zeros(Int, l)
    edg_dsts = zeros(Int, l)
    if args[:bonds]
        edg_lens = zeros(Float64, l)
        for (i, edge) in enumerate(edges(graph))
            edg_srcs[i] = edg_dsts[i + edge_count] = edge.src - 1 # python numbering
            edg_dsts[i] = edg_srcs[i + edge_count] = edge.dst - 1
            edg_lens[i] = edg_lens[i + edge_count] = get_prop(graph, edge, :distance)
        end
        return edg_srcs, edg_dsts, edg_lens
    elseif args[:vspn]
        edg_type = zeros(Int, l)
        edg_wght = zeros(Float64, l)
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


function vspn_feature_matrix(g::MetaGraph, xtal::Crystal, element_to_int::Dict)::SparseMatrixCSC
    nb_atom_types = length(element_to_int)
    embedding_length = nb_atom_types + 1
    X = zeros(Int, nv(g), embedding_length)
    for i in 1:nv(g)
        if get_prop(g, i, :type) == :A # atom
            X[i, element_to_int[xtal.atoms.species[i]]] = 1
            #X[i, nb_atom_types + 1 + degree(xtal.bonds)[i]] = 1
        else # voro pt
            X[i, nb_atom_types + 1] = 1
        end
    end
    return sparse(X)
end


function write_data(xtal::Crystal, name::String, element_to_int::Dict{Symbol,Int}, graphs_path::String, args::Dict{Symbol,Any}, config::Union{Nothing,VSPNConfig}=nothing)
    X_name = chop(name, tail=4)
    if args[:vspn]
        cached("vspn/$X_name.jld2") do
            g = vspn_graph(xtal, config)
            return g, vspn_feature_matrix(g, xtal, element_to_int)
        end
        return
    end
    X = node_feature_matrix(xtal, element_to_int)
	# bond graph
    if args[:bonds]
        A, B, D = edge_vectors(xtal.bonds, args)
        npzwrite(joinpath(graphs_path, X_name * "_node_features.npy"), X)
        npzwrite(joinpath(graphs_path, X_name * "_edges_src.npy"), A)
        npzwrite(joinpath(graphs_path, X_name * "_edges_dst.npy"), B)
        npzwrite(joinpath(graphs_path, X_name * "_euc.npy"), D)
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


function process_example(xtal_name::String, element_to_int::Dict{Symbol,Int}, bonded_xtals_cache::String, graphs_path::String, args::Dict{Symbol,Any}, config::Union{Nothing,VSPNConfig}=nothing)
    @load joinpath(bonded_xtals_cache, xtal_name) obj
    xtal, _ = obj
    write_data(xtal, xtal_name, element_to_int, graphs_path, args, config)
end


function process_examples(good_xtals::Vector{String}, element_to_int::Dict{Symbol,Int}, args::Dict{Symbol,Any})
    if args[:vspn]
        config = VSPNConfig(args[:probe], args[:forcefield])
        @showprogress "Processing examples:" pmap(good_xtals) do x
            process_example(x, element_to_int, rc[:cache][:bonded_xtals], rc[:paths][:graphs], args, config)
        end
    else
        @showprogress "Processing examples:" pmap(good_xtals) do x
            process_example(x, element_to_int, rc[:cache][:bonded_xtals], rc[:paths][:graphs], args)
        end
    end
end
