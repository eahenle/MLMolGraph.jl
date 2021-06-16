Process crystals in `data/crystals` to machine-learning inputs for graph neural networks by:

```
julia --project [-p numprocs] process_crystals.jl [OPTS]
```

Set `-p numprocs` with `numprocs` = desired number of processes for distributed processing.

Processing steps are cached in case of program interruption.  Use `-h` or `--help` to see command line options.
