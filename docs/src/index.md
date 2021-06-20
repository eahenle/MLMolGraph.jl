# MLMolGraph.jl

`MLMolGraph.jl` is a package for converting chemical structures to model inputs for graph neural networks.

Getting started is easy!

1. Install Julia 1.6.0 or higher.

2. In the Julia REPL, hit `]` to enter `Pkg` mode, and install the library:
``` Pkg> add MLMolGraph```

3. In the REPL or a script, import the namespace:
```julia using MLMolGraph```

Once that's done, you're ready to go.  If your use case is covered by the included [processing script](process_crystals), you can even skip step 3.  More about that on the following pages.  If you create your own structure processing pipeline using `MLMolGraph`, please consider making a [pull request](https://github.com/SimonEnsemble/MLMolGraph.jl) to add it to the script!
