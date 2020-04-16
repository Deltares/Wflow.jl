module Wflow

using Dates
using LightGraphs
using NCDatasets
using StaticArrays
using Pkg.TOML

include("config.jl")
include("kinematic_wave.jl")
include("sbm.jl")

end # module
