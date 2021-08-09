module Wflow

using Dates
using TOML
using LightGraphs
using NCDatasets
using StaticArrays
using Statistics
using UnPack
using Random
using BasicModelInterface
using FieldMetadata
using Parameters
using DelimitedFiles
using ProgressLogging
using CFTime
using Base.Threads
using Glob
using NCDatasets: MFDataset

@metadata get_units "mm Δt-1"

const BMI = BasicModelInterface
const Float = Float64
const CFDataset = Union{NCDataset, MFDataset}

mutable struct Clock{T}
    time::T
    iteration::Int
    Δt::Second
end

function Clock(config)
    # this constructor is used by reset_clock!, since if the Clock has already
    # been constructed before, the config is complete
    calendar = get(config, "calendar", "standard")::String
    starttime = cftime(config.starttime, calendar)
    Δt = Second(config.timestepsecs)
    Clock(starttime, 1, Δt)
end

function Clock(config, reader)
    nctimes = reader.dataset["time"][:]
    # if the config file does not have a start or endtime, folow the NetCDF times
    # and add them to the config
    # if the timestep is not given, use the difference between NetCDF time 1 and 2
    starttime = get(config, "starttime", nothing)
    if starttime === nothing
        starttime = first(nctimes)
        config.starttime = starttime
    end
    endtime = get(config, "endtime", nothing)
    if endtime === nothing
        endtime = last(nctimes)
        config.endtime = endtime
    end
    timestepsecs = get(config, "timestepsecs", nothing)
    if timestepsecs === nothing
        timestepsecs = Dates.value(Second(nctimes[2] - nctimes[1]))
        config.timestepsecs = timestepsecs
    end

    calendar = get(config, "calendar", "standard")::String
    starttime = cftime(config.starttime, calendar)
    Δt = Second(timestepsecs)

    fews_run = get(config, "fews_run", false)::Bool
    if fews_run
        starttime = starttime + Δt
    end

    Clock(starttime, 1, Δt)
end

include("io.jl")

"""
    Model{N,L,V,R,W}

Composite type that represents all different aspects of a Wflow Model, such as the
network, parameters, clock, configuration and input and output.
"""
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

include("horizontal_process.jl")
include("hbv.jl")
include("sbm.jl")
include("sediment.jl")
include("reservoir_lake.jl")
include("hbv_model.jl")
include("sbm_model.jl")
include("sediment_model.jl")
include("flow.jl")
include("vertical_process.jl")
include("groundwater/connectivity.jl")
include("groundwater/aquifer.jl")
include("groundwater/boundary_conditions.jl")
include("sbm_gwf_model.jl")
include("utils.jl")
include("bmi.jl")
include("subdomains.jl")

"""
    run(tomlpath::String)
    run(config::Config)
    run(model::Model)
    run()

Run an entire simulation starting either from a path to a TOML settings file,
a prepared `Config` object, or an initialized `Model` object. This allows more flexibility
if you want to for example modify a `Config` before initializing the `Model`.

The 0 argument version expects ARGS to contain a single entry, pointing to the TOML path.
This makes it easier to start a run from the command line without having to escape quotes:

    julia -e "using Wflow; Wflow.run()" "path/to/config.toml"
"""
function run(tomlpath)
    config = Config(tomlpath)
    run(config)
end

function run(config::Config)
    modeltype = config.model.type

    model = if modeltype == "sbm"
        initialize_sbm_model(config)
    elseif modeltype == "sbm_gwf"
        initialize_sbm_gwf_model(config)
    elseif modeltype == "hbv"
        initialize_hbv_model(config)
    elseif modeltype == "sediment"
        initialize_sediment_model(config)
    else
        error("unknown model type")
    end

    run(model)
end

function run(model::Model; close_files = true)
    @unpack network, config, writer, clock = model

    # in the case of sbm_gwf it's currently a bit hard to use dispatch
    model_type = config.model.type::String
    update_func = model_type == "sbm_gwf" ? update_sbm_gwf : update

    # determine timesteps to run
    calendar = get(config, "calendar", "standard")::String
    starttime = clock.time
    Δt = clock.Δt
    endtime = cftime(config.endtime, calendar)
    times = range(starttime, endtime, step = Δt)

    @info "Run information" model_type starttime Δt endtime nthreads()
    @progress for (i, time) in enumerate(times)
        @debug "Starting timestep" time timestep = i
        model = update_func(model)
    end

    # write output state NetCDF
    # undo the clock advance at the end of the last iteration, since there won't
    # be a next step, and then the output state falls on the correct time
    rewind!(clock)
    write_netcdf_timestep(model, writer.state_dataset, writer.state_parameters)

    reset_clock!(model.clock, config)

    # option to support running function twice without re-initializing
    # and thus opening the NetCDF files
    if close_files
        Wflow.close_files(model, delete_output = false)
    end
    return model
end

function run()
    usage = "Usage: julia -e 'using Wflow; Wflow.run()' 'path/to/config.toml'"
    n = length(ARGS)
    if n != 1
        throw(ArgumentError(usage))
    end
    toml_path = only(ARGS)
    if !isfile(toml_path)
        throw(ArgumentError("File not found: $(toml_path)\n" * usage))
    end
    Wflow.run(toml_path)
end

end # module
