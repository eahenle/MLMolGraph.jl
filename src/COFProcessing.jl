module COFProcessing

using CSV, DataFrames, Distributed, JLD2, LightGraphs, MetaGraphs, NPZ, Reexport

@reexport using Xtals


function __init__()
    rc[:cache] = Dict()
    setup_cache()
    rc[:paths][:graphs] = joinpath(rc[:paths][:data], "Graphs")
    if !isdirpath(rc[:paths][:graphs])
        mkpath(rc[:paths][:graphs])
    end
end


function setup_cache()
    rc[:cache][:primitive] = joinpath(rc[:paths][:data], "cache", "primitive")
    rc[:cache][:bonded_xtals] = joinpath(rc[:paths][:data], "cache", "bonded_xtals")
    for dir ∈ rc[:cache]
        if !isdirpath(dir[2]) 
            mkpath(dir[2])
        end
    end
end


function clear_cache()
    try
        rm(joinpath(rc[:paths][:data], "cache"), recursive=true)
    catch
    end
    setup_cache()
end


function cached(f, cachefile)
    cachefile = joinpath(rc[:paths][:data], "cache", cachefile)
    if isfile(cachefile)
        @load cachefile obj
    else
        obj = f()
        @save cachefile obj
    end
    return obj
end


function xtals2primitive(xtal_list)
    @sync @distributed for xtal_file ∈ xtal_list
        cached("primitive/$xtal_file") do 
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


function isgood(xtal_file)
    @load joinpath(rc[:cache][:bonded_xtals], xtal_file) obj
    _, good = obj
    return good
end


function bondNclassify(xtal_list)
    @sync @distributed for xtal_file ∈ xtal_list
        cached("bonded_xtals/$xtal_file") do 
            @load joinpath(rc[:cache][:primitive], xtal_file) obj
            xtal = obj
            if infer_bonds!(xtal, true)
                good = true
            else
                good = false
            end
            return (xtal, good)
        end
    end
    return [xtal_file for xtal_file ∈ xtal_list if isgood(xtal_file)]
end


function encode(xtal_list)
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


function read_targets(source, examples, target_symbol)
    df = CSV.read(joinpath(rc[:paths][:data], source), DataFrame, delim=", ")
    select!(df, [:name, target_symbol])
    filter!(r -> "$(r.name).cif" ∈ examples, df)
    return df
end


function node_feature_matrix(xtal, max_valency, element_to_int)
    embedding_length = length(element_to_int) + max_valency
    X = zeros(Int, xtal.atoms.n, embedding_length)
    for (i, atom) in enumerate(xtal.atoms.species)
        X[i, element_to_int[atom]] = 1
        X[i, embedding_length - degree(xtal.bonds)[i] + 1] = 1
    end
    return X
end


function edge_vectors(graph)
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


function bond_angle_matrix(xtal)
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


function write_data(xtal, name, element_to_int, max_valency)
    X = node_feature_matrix(xtal, max_valency, element_to_int)
    X_name = chop(name, tail=4)
	# bond graph
	A, B, D = edge_vectors(xtal.bonds)
    npzwrite(joinpath(rc[:paths][:graphs], X_name * "_node_features.npy"), X)
    npzwrite(joinpath(rc[:paths][:graphs], X_name * "_edges_src.npy"), A)
    npzwrite(joinpath(rc[:paths][:graphs], X_name * "_edges_dst.npy"), B)
    npzwrite(joinpath(rc[:paths][:graphs], X_name * "_euc.npy"), D)
	# bond angles
    angleABC = bond_angle_matrix(xtal)
    npzwrite(joinpath(rc[:paths][:graphs], X_name * "_angles.npy"), angleABC)
end


function process_example(xtal_name, element_to_int, max_valency)
    @load joinpath(rc[:cache][:bonded_xtals], xtal_name) obj
    xtal, _ = obj
    write_data(xtal, xtal_name, element_to_int, max_valency)
end


function process_examples(good_xtals, element_to_int, max_valency)
    els = [element_to_int for _ ∈ 1:length(good_xtals)]
    mvs = [max_valency for _ ∈ 1:length(good_xtals)]
    pmap(process_example, good_xtals, els, mvs)
end


export cached, xtals2primitive, bondNclassify, encode, read_targets, process_examples, clear_cache

end