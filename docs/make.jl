using Documenter
using MLMolGraph

makedocs(
    root = joinpath(dirname(pathof(Xtals)), "..", "docs"),
    modules = [MLMolGraph],
    sitename = "MLMolGraph.jl",
    clean = true,
    pages = [
            "MLMolGraph" => "index.md",
            "Graphs & AI" => "graph_ml.md",
            "Processing" => "process_crystals.md",
            "API Docs" => "api.md"
            ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(repo = "github.com/SimonEnsemble/MLMolGraph.jl.git")
