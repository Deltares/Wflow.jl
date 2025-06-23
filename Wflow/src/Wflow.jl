module Wflow

import BasicModelInterface as BMI

using Accessors: @optic, @reset, PropertyLens
using Base.Threads: nthreads
using CFTime: CFTime, monthday, dayofyear
using CompositionsBase: decompose
using Dates:
    Dates,
    Second,
    Minute,
    Hour,
    Day,
    Month,
    year,
    TimeType,
    DatePeriod,
    TimePeriod,
    Date,
    DateTime,
    now,
    isleapyear,
    datetime2unix,
    canonicalize
using DelimitedFiles: readdlm
using Glob: glob
using Graphs:
    add_edge!,
    add_vertex!,
    DiGraph,
    dst,
    edges,
    Graph,
    Graphs,
    induced_subgraph,
    inneighbors,
    is_cyclic,
    ne,
    nv,
    outneighbors,
    SimpleDiGraph,
    src,
    topological_sort_by_dfs,
    vertices
using LoggingExtras
using NCDatasets: NCDatasets, NCDataset, dimnames, dimsize, nomissing, defDim, defVar
using Parameters: @with_kw
using Polyester: @batch
using ProgressLogging: @progress
using StaticArrays: SVector, pushfirst, setindex
using Statistics: mean, median, quantile!, quantile
using TerminalLoggers
using TOML: TOML
using Adapt
using KernelAbstractions
import AcceleratedKernels as AK

const CFDataset = Union{NCDataset, NCDatasets.MFDataset}
const CFVariable_MF = Union{NCDatasets.CFVariable, NCDatasets.MFCFVariable}
const VERSION =
    VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
const ROUTING_OPTIONS = (("kinematic-wave", "local-inertial"))

const Float = Float32
const Int = Int32
# const Float = Float64
# const Int = Int64

const cpu_backend = KernelAbstractions.CPU()

using AMDGPU

const gpu_backend = AMDGPU.ROCBackend()
const backend = gpu_backend
const BackendArray = ROCArray
# const backend = cpu_backend
# const BackendArray = Array

# AMDGPU.device!(AMDGPU.devices()[2])

# Allows creation of an array on the backend, e.g. `array_from_host(ones(Int32, 10))` or `array_from_host(1:10, Int32)`
#  copied from https://github.com/JuliaGPU/AcceleratedKernels.jl/blob/main/test/runtests.jl#L140 (MIT License)
function array_from_host(h_arr::AbstractArray, dtype = nothing)
    d_arr = KernelAbstractions.zeros(
        backend,
        isnothing(dtype) ? eltype(h_arr) : dtype,
        size(h_arr),
    )
    copyto!(d_arr, h_arr isa Array ? h_arr : Array(h_arr))      # Allow unmaterialised types, e.g. ranges
    d_arr
end

mutable struct Clock{T}
    time::T
    iteration::Int
    dt::Second
end

function Clock(config)
    # this constructor is used by reset_clock!, since if the Clock has already
    # been constructed before, the config is complete
    calendar = get(config.time, "calendar", "standard")::String
    starttime = cftime(config.time.starttime, calendar)
    dt = Second(config.time.timestepsecs)
    return Clock(starttime, Int(0), dt)
end

function Clock(config, reader)
    nctimes = reader.dataset["time"][:]

    # if the timestep is not given, use the difference between netCDF time 1 and 2
    timestepsecs = get(config.time, "timestepsecs", nothing)
    if timestepsecs === nothing
        timestepsecs = Dates.value(Second(nctimes[2] - nctimes[1]))
        config.time.timestepsecs = timestepsecs
    end
    dt = Second(timestepsecs)

    # if the config file does not have a start or endtime, follow the netCDF times
    # and add them to the config
    starttime = get(config.time, "starttime", nothing)
    if starttime === nothing
        starttime = first(nctimes) - dt
        config.time.starttime = starttime
    end
    endtime = get(config.time, "endtime", nothing)
    if endtime === nothing
        endtime = last(nctimes)
        config.time.endtime = endtime
    end

    calendar = get(config.time, "calendar", "standard")::String
    fews_run = get(config, "fews_run__flag", false)::Bool
    if fews_run
        config.time.starttime = starttime + dt
    end
    starttime = cftime(config.time.starttime, calendar)

    return Clock(starttime, Int(0), dt)
end

abstract type AbstractModel{T} end
abstract type AbstractLandModel end

# different model types (used for dispatch)
abstract type AbstractModelType end
struct SbmModel <: AbstractModelType end         # "sbm" type / sbm_model.jl
struct SbmGwfModel <: AbstractModelType end      # "sbm_gwf" type / sbm_gwf_model.jl
struct SedimentModel <: AbstractModelType end    # "sediment" type / sediment_model.jl

include("io.jl")
include("network.jl")
include("routing/routing.jl")
include("domain.jl")

"""
    Model{R <: Routing, L <: AbstractLandModel, T <: AbstractModelType} <:AbstractModel{T}

Composite type that represents all different aspects of a Wflow Model, such as the network,
parameters, clock, configuration and input and output.
"""
struct Model{R <: Routing, L <: AbstractLandModel, T <: AbstractModelType} <:
       AbstractModel{T}
    config::Config                  # all configuration options
    domain::Domain                  # domain connectivity (network) and shared parameters 
    routing::R                      # routing model (horizontal fluxes), moves along network
    land::L                         # land model simulating vertical fluxes, independent of each other
    clock::Clock                    # to keep track of simulation time
    reader::NCReader                # provides the model with dynamic input
    writer::Writer                  # writes model output
    type::T                         # model type
end

"""
   Model(config::Config)::Model

Initialization of a `Model` based on the `config` object (parsed configuration TOML file
with input, model and output settings).
"""
function Model(config::Config)::Model
    model_type = config.model.type

    if model_type ∉ ("sbm", "sbm_gwf", "sediment")
        error("Unknown model type $model_type.")
    end
    @info "Initialize model variables for model type `$model_type`."

    type = if model_type == "sbm"
        SbmModel()
    elseif model_type == "sbm_gwf"
        SbmGwfModel()
    elseif model_type == "sediment"
        SedimentModel()
    end

    return Model(config, type)
end

# prevent a large printout of model components and arrays
Base.show(io::IO, ::AbstractModel{T}) where {T} = print(io, "model of type ", T)

include("forcing.jl")
include("vegetation/parameters.jl")
include("vegetation/rainfall_interception.jl")
include("vegetation/canopy.jl")
include("snow/snow_process.jl")
include("snow/snow.jl")
include("glacier/glacier_process.jl")
include("glacier/glacier.jl")
include("surfacewater/runoff.jl")
include("soil/soil.jl")
include("soil/soil_process.jl")
include("sbm.jl")
include("groundwater/connectivity.jl")
include("groundwater/aquifer.jl")
include("groundwater/boundary_conditions.jl")
include("routing/timestepping.jl")
include("routing/subsurface.jl")
include("routing/reservoir.jl")
include("routing/lake.jl")
include("routing/surface_kinwave.jl")
include("routing/surface_local_inertial.jl")
include("routing/surface_routing.jl")
include("routing/routing_process.jl")
include("demand/water_demand.jl")
include("sbm_model.jl")
include("sediment/erosion/erosion_process.jl")
include("sediment/erosion/rainfall_erosion.jl")
include("sediment/erosion/overland_flow_erosion.jl")
include("sediment/erosion/soil_erosion.jl")
include("sediment/erosion/river_erosion.jl")
include("sediment/sediment_transport/deposition.jl")
include("sediment/sediment_transport/transport_capacity_process.jl")
include("sediment/sediment_transport/transport_capacity.jl")
include("sediment/sediment_transport/overland_flow_transport.jl")
include("sediment/sediment_transport/land_to_river.jl")
include("sediment/sediment_transport/river_transport.jl")
include("erosion.jl")
include("sediment_flux.jl")
include("sediment_model.jl")
include("routing/initialize_routing.jl")
include("sbm_gwf_model.jl")
include("standard_name.jl")
include("utils.jl")
include("bmi.jl")
include("subdomains.jl")
include("logging.jl")
include("states.jl")

"""
    run(tomlpath::AbstractString; silent=false)
    run(config::Config)
    run!(model::Model)
    run()

Run an entire simulation starting either from a path to a TOML settings file, a prepared
`Config` object, or an initialized `Model` object. This allows more flexibility if you want
to for example modify a `Config` before initializing the `Model`. Logging to a file is only
part of the `run(tomlpath::AbstractString)` method. To avoid logging to the terminal, set
the `silent` keyword argument to `true`, or put that in the TOML.

The 0 argument version expects ARGS to contain a single entry, pointing to the TOML path.
This makes it easier to start a run from the command line without having to escape quotes:

    julia -e "using Wflow; Wflow.run()" "path/to/config.toml"
"""
function run(tomlpath::AbstractString; silent = nothing)
    config = Config(tomlpath)
    # if the silent kwarg is not set, check if it is set in the TOML
    if silent === nothing
        silent = get(config.logging, "silent", false)::Bool
    end
    fews_run = get(config, "fews_run", false)::Bool
    logger, logfile = init_logger(config; silent)
    with_logger(logger) do
        @info "Wflow version `v$VERSION`"
        # to catch stacktraces in the log file a try-catch is required
        try
            run(config)
        catch e
            # avoid logging backtrace for the single line FEWS log format
            # that logger also uses SimpleLogger which doesn't result in a good backtrace
            if fews_run
                @error "Wflow simulation failed" exception = e _id = :wflow_run
            else
                @error "Wflow simulation failed" exception = (e, catch_backtrace()) _id =
                    :wflow_run
            end
            rethrow()
        finally
            close(logfile)
        end
    end
    return nothing
end

function run(config::Config)
    model = Model(config)
    load_fixed_forcing!(model)
    run!(model)
    return model
end

function run_timestep!(model::Model; update_func = update!, write_model_output = true)
    advance!(model.clock)
    load_dynamic_input!(model)
    update_func(model)
    if write_model_output
        write_output(model)
    end
    return nothing
end

function run!(model::Model; close_files = true)
    (; config, writer, clock) = model

    model_type = config.model.type::String

    # determine timesteps to run
    calendar = get(config.time, "calendar", "standard")::String
    starttime = clock.time
    dt = clock.dt
    endtime = cftime(config.time.endtime, calendar)
    times = range(starttime + dt, endtime; step = dt)

    @info "Run information" model_type starttime dt endtime nthreads()
    run_timestep!(model)
    runstart_time = now()

    t0 = time()
    @progress for (i, _time) in enumerate(times)
        if i != 1
            @debug "Starting timestep." _time i now()
            run_timestep!(model)
        end

        if i % 100 == 0
            t1 = time()
            @info "Timestep run time: $((t1-t0)/100)"
            t0 = time()
        end
    end
    @info "Simulation duration: $(canonicalize(now() - runstart_time))"

    # write output state netCDF
    if !isnothing(writer.state_nc_path)
        @info "Write output states to netCDF file `$(writer.state_nc_path)`."
    end
    write_netcdf_timestep(model, writer.state_dataset, writer.state_parameters)

    reset_clock!(model.clock, config)

    # option to support running function twice without re-initializing
    # and thus opening the netCDF files
    if close_files
        Wflow.close_files(model; delete_output = false)
    end

    # copy TOML to dir_output, to archive what settings were used
    if haskey(config, "dir_output")
        src = normpath(pathof(config))
        dst = output_path(config, basename(src))
        if src != dst
            @debug "Copying TOML file." src dst
            cp(src, dst; force = true)
        end
    end
    return nothing
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
    return run(toml_path)
end

end # module
