using Documenter
using MLMolGraph

makedocs(
    root = joinpath(dirname(pathof(MLMolGraph)), "..", "docs"),
    modules = [MLMolGraph],
    sitename = "MLMolGraph",
    clean = true,
    pages = [
        "MLMolGraph" => "index.md",
        "Graphs & AI" => "graph_ml.md"
    ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(repo = "github.com/SimonEnsemble/MLMolGraph.git")
