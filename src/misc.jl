# shifts atomic coordinates to center unit cell on origin and return as matrix
function shift_coords(xtal::Crystal)::Matrix{Float64}
    return Matrix(transpose(Cart(xtal.atoms.coords, xtal.box).x .- xtal.box.f_to_c * [0.5, 0.5, 0.5]))
end

function shift_coords(coords::Matrix{Float64}, box::Box)
    return coords .- box.f_to_c * [0.5, 0.5, 0.5]
end


function shift_back(coords, box)
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
function pbc_distance(pt1::Array{Float64}, pt2::Array{Float64}, box)::Float64
    dx = pt1 .- pt2
    sz = [box.a, box.b, box.c]
    for i ∈ 1:3
        if dx[i] > sz[i] * 0.5
            dx[i] -= sz[i]
        elseif dx[i] ≤ -sz[i] * 0.5
            dx[i] += sz[i]
        end
    end
    return norm(dx)
end


# https://mathworld.wolfram.com/Sphere-SphereIntersection.html
function lens(r, R, d)
    π * (R + r - d)^2 * (d^2 + 2*d*r - 3*r^2 + 2*d*R + 6*r*R - 3*R^2) / 12 / d
end