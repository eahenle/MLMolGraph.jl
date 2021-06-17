testfiles = [
    "voronoi.jl",
    "misc.jl",
    "vspn.jl",
    "cache_tools.jl",
    "processing.jl",
    "run_process.jl"
]
LOGGING_LEVEL = :Debug

using Test, Logging
global_logger(ConsoleLogger(stdout, getproperty(Logging, LOGGING_LEVEL)))

using MLMolGraph
MLMolGraph.banner()

for testfile ∈ testfiles
    @info "Running test/$(testfile)"
    include(testfile)
end
