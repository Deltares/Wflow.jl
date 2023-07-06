# ensure both projects are up to date
using Pkg

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

cd("create_app")

Pkg.activate(".")
Pkg.update()

# remove previous build if it exists
rm("wflow_bundle", force=true, recursive=true)

# build the app, takes ~15 minutes
include("../create_app/create_app.jl")
