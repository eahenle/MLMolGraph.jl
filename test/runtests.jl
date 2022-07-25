testfiles = [
    "processing.jl"
    "run_process.jl"
    "argument_parsing.jl"
    "misc.jl"
]
const DEBUG = false #* toggle debug mode

@assert VERSION.major == 1
@assert VERSION.minor > 6

using Test, Logging
global_logger(ConsoleLogger(stdout, getproperty(Logging, DEBUG ? :Debug : :Info)))
@debug "DEBUG mode"

using MLMolGraph
MLMolGraph.banner()

for testfile ∈ testfiles
    @info "Running test/$(testfile)"
    @time include(testfile)
end

temp_dir = joinpath(pwd(), ".temp")
if isdir(temp_dir)
    rm(temp_dir, recursive=true)
end

@info "Done!"
