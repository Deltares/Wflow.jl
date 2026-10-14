# Usage: julia --project=build/create_binaries build/trim_poc/build.jl <project_dir> <outname>
using JuliaC

project_dir = abspath(ARGS[1])
outname = ARGS[2]
output_dir = abspath(joinpath(@__DIR__, "out", outname))
rm(output_dir; force = true, recursive = true)

ld_flags = String[]
# Use Julia's bundled lld, see build/create_binaries/create_app.jl
if Sys.islinux()
    bundled_lld = normpath(Sys.BINDIR, "..", "libexec", "julia", "lld")
    lld_dir = mktempdir()
    symlink(bundled_lld, normpath(lld_dir, "ld"))
    ld_flags = ["-B$lld_dir"]
end

image_recipe = ImageRecipe(;
    output_type = "--output-exe",
    file = project_dir,
    trim_mode = "safe",
    verbose = true,
)
link_recipe = LinkRecipe(;
    image_recipe,
    outname = joinpath(output_dir, outname),
    rpath = "@bundle",
    ld_flags,
)
bundle_recipe = BundleRecipe(; link_recipe, output_dir)

compile_products(image_recipe)
link_products(link_recipe)
bundle_products(bundle_recipe)
