module Wflow

using Dates
using LightGraphs
using NCDatasets
using StaticArrays
using Statistics
using Setfield: setproperties
using FillArrays
using TypedTables
using UnPack
using CSV
using Random
using CSV.DataFrames

mutable struct Clock
    time::DateTime
    iteration::Int
    Δt::Second
end

struct Model{N,L,V,R,W}
    network::N  # connectivity information, directed graph
    lateral::L  # lateral model that holds lateral state, moves along network
    vertical::V  # vertical model that holds vertical state, independent of each other
    clock::Clock  # to keep track of simulation time
    reader::R  # provides the model with dynamic input
    writer::W  # writes model output
end

include("toml.jl")
include("io.jl")
include("horizontal_process.jl")
include("sbm.jl")
include("reservoir_lake.jl")
include("sbm_model.jl")
include("flow.jl")
include("vertical_process.jl")
include("utils.jl")

end # module
