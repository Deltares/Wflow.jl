using Documenter, Wflow

makedocs(;
    modules=[Wflow],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/Deltares/Wflow.jl/blob/{commit}{path}#L{line}",
    sitename="Wflow.jl",
    authors="Deltares",
    assets=String[],
)

deploydocs(;
    repo="github.com/Deltares/Wflow.jl",
)
