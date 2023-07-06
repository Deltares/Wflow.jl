using PackageCompiler
using TOML
using Markdown

# change directory to this script's location
cd(@__DIR__)

create_app(
    "..",
    "wflow_bundle";
    # map from binary name to julia function name
    executables = ["wflow_cli" => "julia_main"],
    precompile_execution_file = "precompile.jl",
    filter_stdlibs = false  # safer, makes only a tiny difference
)

"""
Add the following metadata files to the newly created build:

- Build.toml
- Project.toml
- Manifest.toml
- README.md
- LICENSE
"""
function add_metadata()
    # save some environment variables in a Build.toml file for debugging purposes
    vars = ["BUILD_NUMBER", "BUILD_VCS_NUMBER"]
    dict = Dict(var => ENV[var] for var in vars if haskey(ENV, var))
    open("wflow_bundle/share/julia/Build.toml", "w") do io
        TOML.print(io, dict)
    end

    # a stripped Project.toml is already added in the same location by PackageCompiler
    # however it is better to copy the original, since it includes the version and compat
    cp("../Project.toml", "wflow_bundle/share/julia/Project.toml", force = true)
    # the Manifest.toml always gives the exact version of Wflow that was built
    cp("../Manifest.toml", "wflow_bundle/share/julia/Manifest.toml", force = true)

    # put the LICENSE in the top level directory
    cp("../LICENSE", "wflow_bundle/LICENSE", force = true)

    # the README.md is also included to the top level directory of the build
    # for a build only the part until the first H2 header is relevant
    # we can use the Markdown stdlib to write only the first part
    # if no H2 header is found, copy the entire README.md
    md = Markdown.parse_file("../README.md")
    h2_begin = findfirst(x -> isa(x, Markdown.Header{2}), md.content)
    if (h2_begin === nothing) || (h2_begin == 1)
        h1_end = length(md.content)
    else
        h1_end = h2_begin - 1
    end
    open("wflow_bundle/README.md", "w") do io
        relevant_part = Markdown.MD(md[1:h1_end])
        print(io, relevant_part)

        # since the exact Wflow version may be hard to find in the Manifest.toml file
        # we can also extract that information, and add it to the README.md
        manifest = TOML.parsefile("../Manifest.toml")
        if !haskey(manifest, "manifest_format")
            error("Manifest.toml is in the old format, run Pkg.upgrade_manifest()")
        end
        julia_version = manifest["julia_version"]
        wflow_entry = only(manifest["deps"]["Wflow"])
        version = wflow_entry["version"]
        # note that that this is not a git commit hash but a git tree hash
        # finding the associated commit is possible, see
        # https://arbitrary-but-fixed.net/git/julia/2021/03/18/git-tree-sha1-to-commit-sha1.html
        # https://discourse.julialang.org/t/git-sha-of-a-project-using-pkg/27211
        tree = wflow_entry["git-tree-sha1"]
        version_info = """

        ## Version

        This wflow_cli build uses the Wflow.jl version mentioned below.

        ```toml
        version = "$version"
        git-tree-sha1 = "$tree"
        julia_version = "$julia_version"
        ```"""
        println(io, version_info)
    end
end

add_metadata()
