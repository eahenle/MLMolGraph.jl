import Pkg

# tries to load a Python dependency; on failure, adds the dependency via Conda
function check_add_dep(pkg; channel="", altname=nothing)
    try
        @info "Checking dependency: $pkg"
        pyimport(isnothing(altname) ? pkg : altname)
    catch
        @info "Installing $pkg..."
        Conda.add(pkg, channel=channel)
        pyimport(isnothing(altname) ? pkg : altname)
        @info "$pkg install verified."
    end
end


@info "Setting up Python environment..."
ENV["PYTHON"]=""
Pkg.add("PyCall")
Pkg.build("PyCall")
using PyCall
pyimport("sys")
@info "PyCall verified."

# check for Conda; if not found, install it
try
    @info "Checking Conda."
    using Conda
    @info "Conda verified."
catch
    @info "Installing Conda..."
    Pkg.add("Conda")
    using Conda
    @info "Conda verified."
end

# check for MLMolGraph; if not found, add it
try
    @info "Checking MLMolGraph."
    using MLMolGraph
    @info "MLMolGraph verified."
catch # if not, install it
    @info "Installing MLMolGraph..."
    Pkg.add("MLMolGraph")
    Pkg.build("MLMolGraph")
    using MLMolGraph
    @info "MLMolGraph install verified."
end

# check the deps, add if missing
check_add_dep("scipy")
check_add_dep("pymatgen", channel="conda-forge")
check_add_dep("pytorch", channel="pytorch", altname="torch")

@info "Setup complete!"
