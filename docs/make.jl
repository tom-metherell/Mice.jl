using Documenter, Mice

makedocs(
    sitename = "Mice.jl",
    modules = [Mice],
    pages = ["index.md",
             "imputation.md",
             "analysis.md",
             "pooling.md",
             "benchmarks.md"]
)