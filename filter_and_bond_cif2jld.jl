using Distributed

@everywhere begin
	using Pkg
	Pkg.activate(".")
	push!(LOAD_PATH, joinpath(pwd(), "src"))
end

@everywhere using Xtals, JLD2

for dir ∈ ["cache", "exclude"]
	path = joinpath(Xtals.PATH_TO_DATA, dir)
	if !isdir(path)
		mkdir(path)
	end
end

@everywhere begin
	xtals = readdir(Xtals.PATH_TO_CRYSTALS)

	function exclude(xtal_file)
		@info "Rejecting $xtal_file"
		mv(joinpath(Xtals.PATH_TO_CRYSTALS, xtal_file), 
			joinpath(Xtals.PATH_TO_DATA, "exclude", xtal_file))
	end

	function filter_and_cache(xtal_file)
		cachefile = joinpath(Xtals.PATH_TO_DATA, "cache", "$xtal_file.jld2")
		try
			if !isfile(cachefile)
				xtal = Crystal(xtal_file, remove_duplicates=true)
				strip_numbers_from_atom_labels!(xtal)
				if infer_bonds!(xtal, true)
					@info "Caching $xtal_file"
					@save cachefile xtal
				else
					exclude(xtal_file)
				end
			else
				@info "Skipping $xtal_file (already cached)"
			end
		catch exception
			@warn "Error processing $xtal_file" exception
			exclude(xtal_file)
		end
	end
end

pmap(filter_and_cache, xtals)
