# Why Graphs?

For many applications, it is advantageous to represent chemical structures as graphs, where the atoms are the graph nodes and the bonds are the graph edges.
This is particularly true in machine learning, where neural networks based on the topology of the chemical graph can be trained to provide on-demand prediction of material properties which would normally require long calculations on powerful computers.

`MLMolGraph` uses the [`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl) package to read chemical structure data and build the bonding graph, and then exports graph information as a series of vectors, formatted on disk for easy application with `torch` or `Flux`.  Training and prediction inputs can readily be produced using the  processing script described on the next page.
