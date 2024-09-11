module Wflow

using Dates
using TOML
using Graphs
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
using LoggingExtras
using TerminalLoggers
using CFTime
using Base.Threads
using Glob
using Polyester
using LoopVectorization
using IfElse

@metadata get_units "mm dt-1" String
@metadata grid_loc "node" String #BMI grid location

const BMI = BasicModelInterface
const Float = Float64
const CFDataset = Union{NCDataset, NCDatasets.MFDataset}
const CFVariable_MF = Union{NCDatasets.CFVariable, NCDatasets.MFCFVariable}
const version =
    VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

mutable struct Clock{T}
    time::T
    iteration::Int
    dt::Second
end

function Clock(config)
    # this constructor is used by reset_clock!, since if the Clock has already
    # been constructed before, the config is complete
    calendar = get(config, "calendar", "standard")::String
    starttime = cftime(config.starttime, calendar)
    dt = Second(config.timestepsecs)
    return Clock(starttime, 0, dt)
end

function Clock(config, reader)
    nctimes = reader.dataset["time"][:]

    # if the timestep is not given, use the difference between netCDF time 1 and 2
    timestepsecs = get(config, "timestepsecs", nothing)
    if timestepsecs === nothing
        timestepsecs = Dates.value(Second(nctimes[2] - nctimes[1]))
        config.timestepsecs = timestepsecs
    end
    dt = Second(timestepsecs)

    # if the config file does not have a start or endtime, follow the netCDF times
    # and add them to the config
    starttime = get(config, "starttime", nothing)
    if starttime === nothing
        starttime = first(nctimes) - dt
        config.starttime = starttime
    end
    endtime = get(config, "endtime", nothing)
    if endtime === nothing
        endtime = last(nctimes)
        config.endtime = endtime
    end

    calendar = get(config, "calendar", "standard")::String
    fews_run = get(config, "fews_run", false)::Bool
    if fews_run
        config.starttime = starttime + dt
    end
    starttime = cftime(config.starttime, calendar)

    return Clock(starttime, 0, dt)
end

include("io.jl")

"""
    Model{N,L,V,R,W}

Composite type that represents all different aspects of a Wflow Model, such as the
network, parameters, clock, configuration and input and output.
"""
struct Model{N, L, V, R, W, T}
    config::Config  # all configuration options
    network::N  # connectivity information, directed graph
    lateral::L  # lateral model that holds lateral state, moves along network
    vertical::V  # vertical model that holds vertical state, independent of each other
    clock::Clock  # to keep track of simulation time
    reader::R  # provides the model with dynamic input
    writer::W  # writes model output
    type::T # model type
end

# different model types (used for dispatch)
struct SbmModel end         # "sbm" type / sbm_model.jl
struct SbmGwfModel end      # "sbm_gwf" type / sbm_gwf_model.jl
struct SedimentModel end    # "sediment" type / sediment_model.jl

# prevent a large printout of model components and arrays
Base.show(io::IO, m::Model) = print(io, "model of type ", typeof(m))

include("forcing.jl")
include("horizontal_process.jl")
include("flow.jl")
include("vegetation/rainfall_interception.jl")
include("vegetation/canopy.jl")
include("snow/snow_process.jl")
include("snow/snow.jl")
include("glacier/glacier_process.jl")
include("glacier/glacier.jl")
include("bucket/bucket.jl")
include("bucket/bucket_process.jl")
include("sbm.jl")
include("demand/water_demand.jl")
include("sediment.jl")
include("reservoir_lake.jl")
include("sbm_model.jl")
include("sediment_model.jl")
include("groundwater/connectivity.jl")
include("groundwater/aquifer.jl")
include("groundwater/boundary_conditions.jl")
include("sbm_gwf_model.jl")
include("utils.jl")
include("bmi.jl")
include("subdomains.jl")
include("logging.jl")
include("states.jl")

"""
    run(tomlpath::AbstractString; silent=false)
    run(config::Config)
    run(model::Model)
    run()

Run an entire simulation starting either from a path to a TOML settings file,
a prepared `Config` object, or an initialized `Model` object. This allows more flexibility
if you want to for example modify a `Config` before initializing the `Model`. Logging to a
file is only part of the `run(tomlpath::AbstractString)` method. To avoid logging to the
terminal, set the `silent` keyword argument to `true`, or put that in the TOML.

The 0 argument version expects ARGS to contain a single entry, pointing to the TOML path.
This makes it easier to start a run from the command line without having to escape quotes:

    julia -e "using Wflow; Wflow.run()" "path/to/config.toml"
"""
function run(tomlpath::AbstractString; silent = nothing)
    config = Config(tomlpath)
    # if the silent kwarg is not set, check if it is set in the TOML
    if silent === nothing
        silent = get(config, "silent", false)::Bool
    end
    fews_run = get(config, "fews_run", false)::Bool
    logger, logfile = init_logger(config; silent)
    with_logger(logger) do
        @info "Wflow version `v$version`"
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
end

function run(config::Config)
    modeltype = config.model.type

    model = if modeltype == "sbm"
        initialize_sbm_model(config)
    elseif modeltype == "sbm_gwf"
        initialize_sbm_gwf_model(config)
    elseif modeltype == "sediment"
        initialize_sediment_model(config)
    else
        error("unknown model type")
    end
    load_fixed_forcing(model)
    return run(model)
end

function run_timestep(model::Model; update_func = update, write_model_output = true)
    advance!(model.clock)
    load_dynamic_input!(model)
    model = update_func(model)
    if write_model_output
        write_output(model)
    end
    return model
end

function run(model::Model; close_files = true)
    @unpack network, config, writer, clock = model

    model_type = config.model.type::String

    # determine timesteps to run
    calendar = get(config, "calendar", "standard")::String
    @warn string(
        "The definition of `starttime` has changed (equal to model state time).\n Please",
        " update your settings TOML file by subtracting one model timestep dt from the",
        " `starttime`, if it was used with a Wflow version up to v0.6.3.",
    )
    starttime = clock.time
    dt = clock.dt
    endtime = cftime(config.endtime, calendar)
    times = range(starttime + dt, endtime; step = dt)

    @info "Run information" model_type starttime dt endtime nthreads()
    runstart_time = now()
    @progress for (i, time) in enumerate(times)
        @debug "Starting timestep." time i now()
        model = run_timestep(model)
    end
    @info "Simulation duration: $(canonicalize(now() - runstart_time))"

    # write output state netCDF
    if haskey(config, "state") && haskey(config.state, "path_output")
        @info "Write output states to netCDF file `$(model.writer.state_nc_path)`."
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
    return run(toml_path)
end

end # module
