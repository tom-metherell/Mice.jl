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
             "benchmarks.md",
             "issues.md",
             "whatsnext.md",
             "acknowledgements.md",
             "references.md"],
    plugins = [bib]
)

deploydocs(
    repo = "https://github.com/tom-metherell/Mice.jl.git",
    versions = "#.#.#"
)