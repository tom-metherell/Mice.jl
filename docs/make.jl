using Documenter, DocumenterCitations, Mice

bib = CitationBibliography("docs/references.bib")

makedocs(
    sitename = "Mice.jl",
    modules = [Mice],
    pages = ["index.md",
             "gettingstarted.md",
             "wrangling.md",
             "imputation.md",
             "analysis.md",
             "pooling.md",
             "multithreading.md",
             "benchmarks.md",
             "issues.md",
             "whatsnext.md",
             "acknowledgements.md",
             "references.md"],
    plugins = [bib]
)

deploydocs(
    repo = "github.com/tom-metherell/Mice.jl.git",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", "v0.0" => "v0.0.0", devurl => devurl]
)
