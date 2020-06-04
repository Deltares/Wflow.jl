using Wflow
using Documenter

makedocs(;
    modules = [Wflow],
    authors = "Deltares and contributors",
    repo = "https://github.com/Deltares/Wflow.jl/blob/{commit}{path}#L{line}",
    sitename = "Wflow.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://Deltares.github.io/Wflow.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/Deltares/Wflow.jl")
