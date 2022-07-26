"""
    dataset = consolidate_data(xtal_names, element_to_int, df, args)
"""
function consolidate_data(good_xtals::Vector{String}, element_to_int::Dict{Symbol, Int}, df::DataFrame, args::Dict{Symbol, Any}, temp_dir::String)::Dict{Symbol, Any}
    # dictionary for collecting data set
    dataset = Dict{Symbol, Any}()
    
    # record args in dictionary
    dataset[:args] = args 
    # collect names into dictionary
    dataset[:names] = good_xtals 
    # collect targets df into dictionary
    dataset[:targets] = df 
    # collect encoding dictionary into dataset (and inverse mapping)
    dataset[:element_to_encoding] = element_to_int
    dataset[:encoding_to_element] = Dict([value => key for (key, value) in element_to_int])

	# collect bond angle vectors into dictionary
    dataset[:bond_angles] = args[:angles] ? (@showprogress "Collecting bond angles " pmap(xtal_name -> load_data(xtal_name, temp_dir)[:angles], good_xtals)) : nothing

    # collect bond distances into dictionary
    #dataset[:bond_distances] = @showprogress "Collecting bond distances " pmap(xtal_name -> load_data(xtal_name, temp_dir)[:bond_edges][3][:], good_xtals)
	
    # collect atom node feature matrices into dictionary
    dataset[:atom_features] = @showprogress "Collecting atom features " pmap(xtal_name -> load_data(xtal_name, temp_dir)[:atom_features], good_xtals)

    ##! distribute process
	dataset[:graph_edges] = @showprogress "Collecting graph edge matrices " [collect_graph_edge_matrix(xtal_name, :bond_graph, temp_dir) for xtal_name in good_xtals]

    dataset[:pvec] = args[:pvec] ? error("##! TODO: pvec implementation") : nothing

    dataset[:elemental_fractions] = @showprogress "Collecting elemental fractions " pmap(xtal_name -> load_data(xtal_name, temp_dir)[:elemental_fractions], good_xtals)
    dataset[:elemental_densities] = @showprogress "Collecting elemental densities " pmap(xtal_name -> load_data(xtal_name, temp_dir)[:elemental_densities], good_xtals)
   
    return dataset
end


function reprocess((gem, atom_x))
    if isnothing(gem)
        return nothing, nothing, nothing
    end

    torch = rc[:torch]
    numpy = rc[:numpy]

    nb_atoms = atom_x.size(0)

    # pull out the src/dst lists
    atom_edge_src = [gem[1][i] for i in eachindex(gem[1]) if gem[3][i] == 0]
    atom_edge_dst = [gem[2][i] for i in eachindex(gem[1]) if gem[3][i] == 0]

    # convert directed edges to undirected
    atom_edge_src, atom_edge_dst = vcat(atom_edge_src, atom_edge_dst), vcat(atom_edge_dst, atom_edge_src)

    # convert src/dst arrays to edge_index tensor
    atom_edge_index = torch.tensor(numpy.array([atom_edge_src, atom_edge_dst]), dtype=torch.long)

    return atom_edge_index
end


function export_data(dataset::Dict{Symbol, Any}, args::Dict)
    torch  = rc[:torch]
    numpy  = rc[:numpy]
    pickle = rc[:pickle]

    xtal_name = dataset[:names]
	
	atom_x = !isnothing(dataset[:atom_features]) ? [torch.tensor(atom_node_fts, dtype=torch.float) for atom_node_fts in dataset[:atom_features]] : nothing
	
	graph_edge_matrix = !isnothing(dataset[:graph_edges]) ? [[gem[:, 1] .- 1, gem[:, 2] .- 1, gem[:, 3] .- 1] for gem in dataset[:graph_edges]] : nothing

    atom_edge_index = @showprogress "Finalizing edge index list " pmap(reprocess, zip(graph_edge_matrix, atom_x))
	
	x_p = !isnothing(dataset[:pvec]) ? [torch.tensor(numpy.array(pvec), dtype=torch.float) for pvec in dataset[:pvec]] : nothing
	
	y = !isnothing(dataset[:targets]) ? torch.tensor(numpy.array(dataset[:targets][:, dataset[:args][:target]]), dtype=torch.float) : nothing

	atom_to_int = dataset[:element_to_encoding]

    elemental_densities = dataset[:elemental_densities]
    elemental_fractions = dataset[:elemental_fractions]

    data_dict = Dict([
        "xtal_name"         => xtal_name
        "atom_x"            => atom_x
        "x_p"               => x_p
        "y"                 => y
		"atom_to_int"       => atom_to_int
		"atom_edge_index"   => atom_edge_index
        "elemental_densities" => elemental_densities
        "elemental_fractions" => elemental_fractions
    ])

    # write dataset dictionary to python-readable binary
    PyCall.open(args[:dataset_output], "w") do file
		pickle.dump(data_dict, file)
	end
end


"""
    data = load_data(xtal_name, temp_dir)
Loads the data dictionary for the named input (usually a crystal identifier) from the temporary data store at temp_dir
"""
function load_data(name::String, temp_dir::String)::Dict
    return load_object(joinpath(temp_dir, name))
end


"""
    save_data(name, data, temp_dir)
Write a data dictionary to the temporary data store at temp_dir with the specified name
"""
function save_data(name::String, data::Dict, temp_dir::String)
    save_object(joinpath(temp_dir, name), data)
end


function collect_graph_edge_matrix(xtal_name::String, graphkey::Symbol, temp_dir::String)::Matrix{Int}
	graph = load_data(xtal_name, temp_dir)[graphkey]
	gem = zeros(Int, (3, ne(graph)))
	edge_label = Dict([
		:AA => 1
	])
	for (i, e) in enumerate(edges(graph))
		gem[:, i] .= src(e), dst(e), :type in keys(props(graph, e)) ? edge_label[get_prop(graph, e, :type)] : edge_label[:AA]
	end

	return Matrix(gem')
end
