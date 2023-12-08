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
             "rcall.md",
             "multithreading.md",
             "benchmarks.md",
             "issues.md",
             "whatsnext.md",
             "acknowledgements.md",
             "references.md"],
    plugins = [bib],
    format = Documenter.HTML(
        analytics = "G-0SJ9WPE2ZH"
    )
)

deploydocs(
    repo = "github.com/tom-metherell/Mice.jl.git",
    versions = ["stable" => "v^", "v#.#", "v0.0" => "v0.0.0", "dev" => "dev"]
)
