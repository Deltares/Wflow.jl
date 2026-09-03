using JuliaC
import PackageCompiler
using PackageCompiler: default_app_cpu_target
using TOML
using Pkg
using LicenseCheck

# change directory to this script's location
cd(@__DIR__)

project_dir = "../wflow_cli"
license_file = "../../LICENSE"
output_dir = "wflow_bundle"
git_repo = "../.."
sbom_file = "../../Wflow.spdx.json"

rm(output_dir; force = true, recursive = true)

image_recipe = ImageRecipe(;
    output_type = "--output-exe",
    file = joinpath(project_dir, "src", "wflow_cli.jl"),
    project = project_dir,
    cpu_target = default_app_cpu_target(),
    verbose = true,
)
link_recipe = LinkRecipe(;
    image_recipe,
    outname = joinpath(output_dir, "wflow_cli"),
    rpath = "@bundle",
)
bundle_recipe = BundleRecipe(; link_recipe, output_dir)

precompile_config_dir = abspath("../../Wflow/test")
withenv("WFLOW_PRECOMPILE_CONFIG_DIR" => precompile_config_dir) do
    compile_products(image_recipe)
end
link_products(link_recipe)
bundle_products(bundle_recipe)

include("add_metadata.jl")
add_metadata(project_dir, license_file, output_dir, git_repo, sbom_file)
