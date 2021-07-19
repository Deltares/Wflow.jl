using Wflow
using Documenter
using Plots
using LaTeXStrings

pages = [
    "Home" => "index.md",
    "Publications" => "publications.md",
    "Quick start" => "quick-start.md",
    "Command line interface" => "cli.md",
    "Structures and config" =>
        ["Structures" => "structure.md", "Config" => "config.md"],
    "Models" => [
        "wflow_sbm" => "model/sbm.md",
        "wflow_hbv" => "model/hbv.md",
        "wflow_sediment" => "model/sediment.md",
    ],
    "Vertical components" => [
        "SBM" => "vertical/sbm.md",
        "HBV" => "vertical/hbv.md",
        "Vertical processes" => "vertical/process.md",
        "Soil Loss" => "vertical/sediment.md",
    ],
    "Lateral components" => [
        "Kinematic wave" => "lateral/kinwave.md",
        "Groundwater flow" => "lateral/gwf.md",
        "Sediment Flux" => "lateral/sediment.md",
    ],
    "Basic Model Interface" => "bmi.md",
    "Run from Delft-FEWS" => "fews.md",
    "Changelog" => "changelog.md",
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
    pages = pages,
)

deploydocs(; repo = "github.com/Deltares/Wflow.jl", push_preview = true)
