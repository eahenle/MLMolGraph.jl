# calculates distance between two points in a periodic system
function pbc_distance(pt1::Array{Float64}, pt2::Array{Float64}, box::Vector{Float64})::Float64
    dx = pt1 .- pt2
    for i ∈ 1:3
        if dx[i] > box[i] * 0.5
            dx[i] -= box[i]
        elseif dx[i] ≤ -box[i] * 0.5
            dx[i] += box[i]
        end
    end
    return norm(dx)
end

pbc_distance(pt1::Array{Float64}, pt2::Array{Float64}, box::Box) = pbc_distance(pt1, pt2, [box.a, box.b, box.c])


function validate_args(args::Dict{Symbol,Any})::Tuple{Bool,String}
    if !ispath(args[:data])
        return false, "Invalid data path"
    end
    if args[:angles] && !args[:bonds]
        return false, "Must specify --bonds when requesting --angles"
    end
    return true, ""
end


function banner()
    FIGlet.render("MLMolGraph", FIGlet.availablefonts()[441])
end


_primitive_xtal_name(xtal_name::String) = split(xtal_name, ".cif")[1] * "_primitive_cell.cif"


function record_args(args::Dict{Symbol, Any})
    open(joinpath(rc[:paths][:data], "data_processing_args.csv"), "w") do file
        for (key, value) in args
            write(file, "$key, $value\n")
        end
    end
end
