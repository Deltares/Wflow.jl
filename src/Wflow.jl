module Wflow

using Dates
using LightGraphs
using NCDatasets
using Pkg.TOML
using StaticArrays
using Statistics
using Setfield: setproperties

include("config.jl")
include("horizontal_process.jl")
include("kinematic_wave.jl")
include("model.jl")
include("sbm.jl")
include("sbm_model.jl")
include("subsurface_flow.jl")
include("vertical_process.jl")
include("utils.jl")
include("model.jl")

end # module
