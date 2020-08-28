module Wflow

using Dates
using LightGraphs
using NCDatasets
using StaticArrays
using Statistics
using Setfield: setproperties
using UnPack
using CSV
using Random
using CSV.DataFrames

mutable struct Clock
    time::DateTime
    iteration::Int
    Δt::Second
end

include("io.jl")

struct Model{N,L,V,R,W}
    config::Config  # all configuration options
    network::N  # connectivity information, directed graph
    lateral::L  # lateral model that holds lateral state, moves along network
    vertical::V  # vertical model that holds vertical state, independent of each other
    clock::Clock  # to keep track of simulation time
    reader::R  # provides the model with dynamic input
    writer::W  # writes model output
end

# prevent a large printout of model components and arrays
Base.show(io::IO, m::Model) = print(io, "model of type ", typeof(m))

include("toml.jl")
include("name.jl")
include("horizontal_process.jl")
include("sbm.jl")
include("reservoir_lake.jl")
include("sbm_model.jl")
include("flow.jl")
include("vertical_process.jl")
include("utils.jl")

"""
    run_simulation(tomlpath::String)
    run_simulation(config::Config)
    run_simulation(model::Model)

Run an entire simulation starting either from a path to a TOML settings file,
a prepared `Config` object, or an initialized `Model` object. This allows more flexibility
if you want to for example modify a `Config` before initializing the `Model`.
"""
function run_simulation(tomlpath)
    config = Wflow.Config(tomlpath)
    run_simulation(config)
end

function run_simulation(config::Config)
    model = Wflow.initialize_sbm_model(config)
    run_simulation(model)
end

function run_simulation(model::Model; close_files = true)
    @unpack network, config = model

    times = config.starttime:Second(config.timestepsecs):config.endtime
    for _ in times
        model = Wflow.update(model)
    end

    reset_clock!(model.clock, config)

    # option to support running function twice without re-initializing
    # and thus opening the NetCDF files
    if close_files
        Wflow.close_files(model, delete_output = false)
    end
    return model
end

end # module
