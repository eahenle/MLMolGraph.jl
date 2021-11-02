# shifts atomic coordinates to center unit cell on origin and return as matrix
function shift_coords(xtal::Crystal)::Matrix{Float64}
    return Matrix(transpose(Cart(xtal.atoms.coords, xtal.box).x .- xtal.box.f_to_c * [0.5, 0.5, 0.5]))
end

function shift_coords(coords::Matrix{Float64}, box::Box)::Matrix{Float64}
    return coords .- box.f_to_c * [0.5, 0.5, 0.5]
end


function shift_back(coords::Vector{Float64}, box::Box)::Vector{Float64}
    return coords .+ box.f_to_c * [0.5, 0.5, 0.5]
end


# converts an Xtals box to a Freud box (unit cell tilt factors)
function convert_box(box::Box)::Array{Float64,1}
    v1 = [box.a, 0, 0]
    v2 = [box.b * cos(box.γ), box.b, 0]
    v3 = [box.c * cos(box.β), box.c * cos(box.α), box.c]
    # math from https://hoomd-blue.readthedocs.io/en/stable/box.html
    Lx = norm(v1)
    a2x = dot(v1, v2) / Lx
    Ly = √(norm(v2)^2 - a2x^2)
    xy = a2x / Ly
    Lz = dot(v3, (cross(v1, v2)) / norm(cross(v1, v2)))
    a3x = dot(v1, v3) / Lx
    xz = a3x / Lz
    yz = (dot(v2, v3) - a2x * a3x) / (Ly * Lz)
    return [box.a, box.b, box.c, xy, xz, yz]
end


# calculates distance between two points in a periodic system
## is this correct for non-orthorhombic cells?
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


# https://mathworld.wolfram.com/Sphere-SphereIntersection.html
function lens(r::Float64, R::Float64, d::Float64)::Float64
    π * (R + r - d)^2 * (d^2 + 2*d*r - 3*r^2 + 2*d*R + 6*r*R - 3*R^2) / 12 / d
end


function validate_args(args::Dict{Symbol,Any})::Tuple{Bool,String}
    if !ispath(args[:data])
        return false, "Invalid data path"
    end
    if args[:angles] && !args[:bonds]
        return false, "Must specify bonds when requesting angles"
    end
    if args[:vspn] && (args[:forcefield] == "" || args[:probe] == "")
        return false, "Must specify probe and forcefield for VSPN"
    end
    if !ispath(args[:cache])
        return false, "Invalid cache path"
    end
    return true, ""
end


function banner()
    FIGlet.render("MLMolGraph", FIGlet.availablefonts()[64])
end
