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

ld_flags = String[]
# Workaround for `lib/julia/libgcc.a: error adding symbols: File format not recognized`
# Use Julia's bundled lld, to avoid old host ld not being able to read the compressed file
# See https://github.com/JuliaLang/julia/pull/61652
if Sys.islinux()
    bundled_lld = normpath(Sys.BINDIR, "..", "libexec", "julia", "lld")
    isfile(bundled_lld) || error("Julia's bundled lld linker was not found at $bundled_lld")
    lld_dir = mktempdir()
    symlink(bundled_lld, normpath(lld_dir, "ld"))
    ld_flags = ["-B$lld_dir"]
end

# JuliaC links the executable directly instead of loading a relocatable sysimage.
image_recipe = ImageRecipe(;
    output_type = "--output-exe",
    file = project_dir,
    cpu_target = default_app_cpu_target(),
    verbose = true,
)
link_recipe = LinkRecipe(;
    image_recipe,
    outname = joinpath(output_dir, "wflow_cli"),
    rpath = "@bundle",
    ld_flags,
)
bundle_recipe = BundleRecipe(; link_recipe, output_dir)

compile_products(image_recipe)
link_products(link_recipe)
bundle_products(bundle_recipe)

include("add_metadata.jl")
add_metadata(project_dir, license_file, output_dir, git_repo, sbom_file)
