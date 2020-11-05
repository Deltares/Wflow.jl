using Wflow
using Documenter

pages = [
    "Home" => "index.md",
    "Quick start" => "quick-start.md",
    "Model structure" => "structure.md",
    "Models" => [
        "wflow_sbm" => "model/sbm.md",
        "wflow_hbv" => "model/hbv.md",
    ],
    "Vertical components" => [
        "SBM" => "vertical/sbm.md",
        "HBV" => "vertical/hbv.md",
        "Vertical processes" => "vertical/process.md",
    ],
    "Lateral components" => [
        "Kinematic wave" => "lateral/kinwave.md",
        "Groundwater flow" => "lateral/gwf.md",
    ],
    "Developer documentation" => "devdocs.md",
]

makedocs(;
    modules = [Wflow],
    authors = "Deltares and contributors",
    repo = "https://github.com/Deltares/Wflow.jl/blob/{commit}{path}#L{line}",
    sitename = "Wflow.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://deltares.github.io/Wflow.jl",
        assets = String[],
    ),
    pages = pages
)

deploydocs(; repo = "github.com/Deltares/Wflow.jl")
