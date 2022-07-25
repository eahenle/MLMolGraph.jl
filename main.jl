DEBUG = false # use to toggle debugging outputs

using Distributed, Random
@everywhere using Logging

@everywhere begin ## TODO remove pre-registration dev hack
    global_logger(ConsoleLogger(stdout, Logging.Error))
    import Pkg
    Pkg.activate(".")
    global_logger(ConsoleLogger(stdout, Logging.Info))
end

global_logger(Logging.ConsoleLogger(stdout, DEBUG ? Logging.Debug : Logging.Info))
@debug "DEBUG mode"

@everywhere using MLMolGraph

MLMolGraph.banner()

# read the command line arguments
@everywhere args = parse_args()

# check args
valid, msg = MLMolGraph.validate_args(args)
if !valid
    error(msg)
end

@info "Program parameters" args...

# set the random seed (seeded randomly by default)
Random.seed!(args[:seed])

# set paths relative to the input data folder
set_paths(args[:data])

# run the data processing
run_process(args)

@info "Done!"
exit(0) ##! why does the end of the script throw some silly error when nothing has failed?
