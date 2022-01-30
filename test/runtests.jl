testfiles = [
    "misc.jl",
    "cache_tools.jl",
    "voronoi.jl",
    "processing.jl", 
    "run_process.jl"
]
LOGGING_LEVEL = :Info

@assert VERSION.major == 1
@assert VERSION.minor > 5

using Test, Logging
global_logger(ConsoleLogger(stdout, getproperty(Logging, LOGGING_LEVEL)))

using MLMolGraph
MLMolGraph.banner()

for testfile ∈ testfiles
    @info "Running test/$(testfile)"
    @time include(testfile)
end

@info "Done!"
