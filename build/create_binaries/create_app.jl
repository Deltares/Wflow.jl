using PackageCompiler
using TOML
using LibGit2
using Pkg

# change directory to this script's location
cd(@__DIR__)

project_dir = "../wflow_cli"
license_file = "../../LICENSE"
output_dir = "wflow_bundle"
git_repo = "../.."

# build the latest release by default, unless WFLOW_REV is set
# WFLOW_REV is an optional git revision. rev can be a branch name or a git commit SHA1
build_rev = haskey(ENV, "WFLOW_REV")

Pkg.activate(".")
if build_rev
    rev = ENV["WFLOW_REV"]
    Pkg.add(url="https://github.com/Deltares/Wflow.jl", rev=rev)
end
Pkg.update()

# confirm Wflow passes tests and download the required data to run precompile.jl
using Wflow
Pkg.test("Wflow")

# cd("create_app")

Pkg.activate(".")
Pkg.update()

# remove previous build if it exists
rm("output_dir", force=true, recursive=true)

create_app(
    project_dir,
    output_dir;
    # map from binary name to julia function name
    executables=["wflow_cli" => "julia_main"],
    precompile_execution_file="precompile.jl",
    filter_stdlibs=false,
    force=true
)

include("add_metadata.jl")
add_metadata(project_dir, license_file, output_dir, git_repo)
