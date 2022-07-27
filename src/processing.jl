using XtalsPyTools

function loadcheck(xtal_list::Vector{String}, tolerance::Float64, temp_dir::String)
    @showprogress "Checking file loadability. " @distributed for xtal_name ∈ xtal_list
        xtal_data = Dict{Symbol, Any}()
        try
            xtal_data[:input] = Crystal(xtal_name * ".cif", remove_duplicates=true, net_charge_tol=tolerance)
            save_data(xtal_name, xtal_data, temp_dir)
        catch exception
            @error xtal_name exception
            save_data(xtal_name, xtal_data, temp_dir)
        end
    end
end


function xtals2primitive(xtal_list::Vector{String}, temp_dir::String)
    @showprogress "Converting to primitive cells. " @distributed for xtal_name ∈ xtal_list
        xtal_data = load_data(xtal_name, temp_dir)
        xtal = xtal_data[:input]
        if isnothing(xtal)
            continue
        else
            xtal_data[:primitive_cell] = primitive_cell(xtal)
            save_data(xtal_name, xtal_data, temp_dir)
        end
    end    
end


function isgood(xtal_file::String, temp_dir)::Bool
    try
        xtal_data = load_data(xtal_file, temp_dir)
        return !isnothing(xtal_data[:bond_graph]) && ne(xtal_data[:bond_graph]) > 0
    catch
        return false
    end
end


function bondNclassify(xtal_list::Vector{String}, temp_dir::String)::Vector{String}
    @showprogress "Inferring bonds. " pmap(xtal_list) do xtal_name
        try
            xtal_data = load_data(xtal_name, temp_dir)
            xtal = xtal_data[:primitive_cell]
            made_bonds = infer_bonds!(xtal, true, calculate_vectors=true)
            if made_bonds
                xtal_data[:bond_graph] = xtal.bonds
            else
                xtal_data[:bond_graph] = nothing
            end
            save_data(xtal_name, xtal_data, temp_dir)
        catch exception
            @error "Exception while classifying input" xtal_name exception
        end
    end
    good = SharedArray{Bool}(length(xtal_list))
    @sync @distributed for i ∈ eachindex(xtal_list)
        good[i] = isgood(xtal_list[i], temp_dir)
    end
    return xtal_list[good]
end


function unique_elements(xtal_file::String, args, temp_dir::String)
    xtal_data = load_data(xtal_file, temp_dir)
    elements = unique(xtal_data[:primitive_cell].atoms.species)
    if args[:verbose]
        @info "$elements unique elements in $(xtal_file)"
    end
    return elements
end


function determine_encoding(xtal_list::Vector{String}, args, temp_dir::String)
    n = length(xtal_list)
    all_elements = [keys(rc[:covalent_radii])...] # list of all elements known to Xtals
    elements = SharedArray{Bool}(length(all_elements)) # elements[i] == true iff keys(all_elements)[i] ∈ some xtal
    if args[:verbose]
        @info "Encoding $n structures' atoms"
        @info "$(length(all_elements)) elements"
    end
    @showprogress "Determining encoding scheme:" @distributed for i ∈ 1:n
        xtal_elements = unique_elements(xtal_list[i], args, temp_dir)
        elements[[findfirst(isequal(element), all_elements) for element in xtal_elements]] .= true
    end
    element_to_int = Dict{Symbol,Int}([element => i for (i, element) ∈ enumerate(all_elements[elements])])
    return element_to_int, length(element_to_int)
end


function read_targets(source::String, examples::Vector{String}, target_symbol::Symbol)::DataFrame
    df = CSV.read(joinpath(rc[:paths][:data], source), DataFrame, delim=",")
    select!(df, [:name, target_symbol])
    filter!(r -> "$(r.name).cif" ∈ examples, df)
    return df
end


# one-hot encoding of atomic species
function node_feature_matrix(xtal::Crystal, element_to_int::Dict{Symbol,Int})
    X = zeros(Int, xtal.atoms.n, length(element_to_int))
    for (i, atom) in enumerate(xtal.atoms.species)
        X[i, element_to_int[atom]] = 1
    end
    return X
end


function edge_vectors(graph::MetaGraph, type::Symbol)
    edge_count = ne(graph)
    l = 2 * edge_count
    edg_srcs = zeros(Int, l)
    edg_dsts = zeros(Int, l)
    edg_prop = zeros(Float64, l)
    for (i, edge) in enumerate(edges(graph))
        edg_srcs[i] = edg_dsts[i + edge_count] = edge.src
        edg_dsts[i] = edg_srcs[i + edge_count] = edge.dst
        edg_prop[i] = edg_prop[i + edge_count] = get_prop(graph, edge, :distance)
    end
    return edg_srcs, edg_dsts, edg_prop
end


function calculate_elemental_fractions(xtal::Crystal, element_to_int::Dict{Symbol,Int})::Vector{Float64}
	fractions = zeros(length(element_to_int))
	ef = empirical_formula(xtal)
    n = sum(values(ef))
    fractions[[element_to_int[x] for x in Symbol.(keys(ef))]] .= [y/n for y in values(ef)]
	return fractions
end


function calculate_elemental_densities(xtal::Crystal, element_to_int::Dict{Symbol,Int})::Vector{Float64}
	densities = zeros(length(element_to_int))
	ef = empirical_formula(xtal)
    densities[[element_to_int[x] for x in Symbol.(keys(ef))]] .= values(ef) ./ sum(values(empirical_formula(xtal))) / xtal.box.Ω
	return densities
end


function write_data(xtal_data::Dict, name::String, element_to_int::Dict{Symbol,Int}, df::DataFrame, args::Dict{Symbol,Any}, temp_dir::String)
	xtal = xtal_data[:primitive_cell]
    # node features
    xtal_data[:atom_features] = node_feature_matrix(xtal, element_to_int)
    # bond graph
    xtal_data[:bond_edges] = args[:bonds] ? edge_vectors(xtal_data[:bond_graph], :bonds) : nothing
    # pvec
    xtal_data[:pvec] = args[:pvec] ? [df[findfirst(n -> n == name, df.name), prop] for prop in args[:pvec_columns]] : nothing
    # element density matrices
    xtal_data[:elemental_fractions] = calculate_elemental_fractions(xtal, element_to_int)
    xtal_data[:elemental_densities] = calculate_elemental_densities(xtal, element_to_int)

    save_data(name, xtal_data, temp_dir)
end


function process_examples(good_xtals::Vector{String}, element_to_int::Dict{Symbol,Int}, df::DataFrame, args::Dict{Symbol,Any}, temp_dir::String)::Vector{Bool}
    @assert length(good_xtals) > 0 "No inputs"
    @showprogress "Processing examples " pmap(good_xtals) do xtal_name
        write_data(load_data(xtal_name, temp_dir), xtal_name, element_to_int, df, args, temp_dir)
    end
    good = trues(length(good_xtals))
    return good
end
