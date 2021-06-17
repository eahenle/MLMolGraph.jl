Process crystals in `data/crystals` to machine-learning inputs for graph neural networks by:

```
julia --project [-p numprocs] process_crystals.jl [OPTS]
```

Set `-p numprocs` with `numprocs` = desired number of processes for distributed processing.

Processing steps are cached in case of program interruption.  Use `-h` or `--help` to see command line options.

**Example uses**

Generate ML inputs using 12 cores for an model based on the bonding network with angle information:

```
julia --project -p 12 process_crystals.jl --bonds --angles
```

Generate ML inputs using 32 cores for a model using VSPN with methane and UFF:

```
julia --project -p 32 process_crystals.jl --vspn -p CH4 -f UFF
```
