testfiles = [
    "misc.jl",
    "cache_tools.jl",
    "voronoi.jl"#,
    ## pymatgen doesn't play nice w/ these inputs...
    #"processing.jl", 
    #"run_process.jl"
]
LOGGING_LEVEL = :Info

using Test, Logging
global_logger(ConsoleLogger(stdout, getproperty(Logging, LOGGING_LEVEL)))

using MLMolGraph
MLMolGraph.banner()

for testfile ∈ testfiles
    @info "Running test/$(testfile)"
    include(testfile)
end
