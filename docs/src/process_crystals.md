# process_crystals.jl

The processing script in the root of the GitHub repository allows the user to process large data sets with nearly zero setup overhead and no custom code.
Simply follow the library installation directions on the landing page (recap: `pkg> add MLMolGraph`) and run the script like so:

`julia [-p numprocs] process_crystals.jl [OPTS]`

`numprocs` is the number of processors to use.  For large data sets, it may be desirable to use the maximum number of available CPUs.
To see the full list of run options, use the `--help` or `-h` flag:

```
$ julia process_crystals.jl --help
usage: process_crystals [-t TARGET] [-d DATA] [-b] [-a] [-v]
                        [-f FORCEFIELD] [-p PROBE] [-x] [-c CACHE]
                        [-h]

reads .cif inputs, converts to primitive cell, infers bonds, and
converts valid bonded structures to ML inputs for GNNs

optional arguments:
  -t, --target TARGET   training target (type: Symbol, default:
                        Symbol("deliverable capacity [v STP/v]"))
  -d, --data DATA       root of data tree. inputs loaded from
                        `crystals` subdirectory. (default:
                        "/home/user/project/data")
  -b, --bonds           write ML inputs for the bonding graph
  -a, --angles          write ML inputs for bonding angles
  -v, --vspn            write ML inputs for Voronoi sphere pore
                        network
  -f, --forcefield FORCEFIELD
                        the forcefield to use for VSPN (default: "")
  -p, --probe PROBE     the molecule to use as the probe in VSPN
                        (default: "")
  -x, --clear_cache     delete the cache
  -c, --cache CACHE     cache location (default:
                        "/home/user/project/data/cache")
  -h, --help            show this help message and exit
```

So, for example, if you wanted to train a model on a data set located at `~/my_data/`, using bonding information including angles and distances, to predict the property `target1`, and had 8 cores available to process the data, you might do like so:

`julia -p 8 process_crystals.jl -d ~/my_data/ -t target1 -b -a`

Of course, many novel and exciting machine learning models require different data to be encoded than what are offered with this program.
The following page details the `MLMolGraph` API to assist in adapting the code for new purposes.
