using Wflow
using Documenter

pages = [
    "Introduction" => [
        "index.md",
        # "intro/use_cases.md",
        "intro/publications.md",
    ],
    "User guide" => [
        "user_guide/intro.md",
        "user_guide/install.md",
        "How to use wflow" => [
            "user_guide/step1_requirements.md",
            "user_guide/step2_settings_file.md",
            "user_guide/step3_input_data.md",
            "user_guide/step4_running.md",
            "user_guide/step5_output.md",
            "user_guide/additional_options.md",
            "user_guide/sample_data.md",
            ],
        "user_guide/model-setup.md",

    ],
    "Model documentation" => [
        "model_docs/intro.md",
        "Vertical concepts" => [
            # "model_docs/vertical/sbm.md",
            "model_docs/vertical/sbm_new.md",
            "model_docs/vertical/hbv.md",
            "model_docs/vertical/sediment.md",
            "model_docs/shared_concepts.md",
        ],
        "Lateral concepts" => [
            "model_docs/lateral/gwf.md",
            "model_docs/lateral/kinwave.md",
            "model_docs/lateral/sediment_flux.md",
            
        ],
        "Model parameters" => [
            "model_docs/params_vertical.md",
            "model_docs/params_lateral.md",
            "model_docs/structures.md",
        ]
    ],       
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
		collapselevel = 2,
    ),
    pages = pages,
)

deploydocs(; repo = "github.com/Deltares/Wflow.jl", push_preview = true)