| **Documentation** | **Build Status** | **Test Coverage** |
|:---:|:---:|:---:|
| [![Docs]()]() | [![Build](https://github.com/eahenle/MLMolGraph.jl/actions/workflows/ci_testing.yml/badge.svg)](https://github.com/eahenle/MLMolGraph.jl/actions/workflows/ci_testing.yml) | [![codecov]()]() [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) |

# Data Processing

Process crystals in `data/crystals` to machine-learning inputs for graph neural networks by:

```
julia --project [-p numprocs] process_crystals.jl [OPTS]
```

Set `-p numprocs` with `numprocs` = desired number of processes for distributed processing.

Use `-h` or `--help` to see command line options.

**Example uses**

Generate ML inputs using 12 cores for a model based on the bonding network with angle information:

```
julia --project -p 12 process_crystals.jl --bonds --angles
```
