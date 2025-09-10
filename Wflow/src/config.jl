#=
Code to support input and output of data and configuration.
Input data can be loaded from netCDF files.
Output data can be written to netCDF or CSV files.
For configuration files we use TOML.
=#

abstract type AbstractConfigSection end

# Don't error on extra fields that aren't hardcoded below
Configurations.ignore_extra(::Type{<:AbstractConfigSection}) = true

# Configurations.jl reads fields of this type as Dict{String, Any}.
# These are subsequently converted to (nested) PropertyDict objects
# for easy accessing of the data.
const MaybePropertyDict =
    Union{Dict{String, Any}, PropertyDict{String, Any, Dict{String, Any}}}

# Pretty printing AbstractConfigurationSection instances
function Base.show(io::IO, c::AbstractConfigSection)
    first = true
    for field in fieldnames(typeof(c))
        f = getfield(c, field)
        if f !== nothing
            first && (first = false; println(io))
            println(io, "\t\t$field\t= $f")
        end
    end
end

# Time related configurations
@option mutable struct TimeSection <: AbstractConfigSection
    const calendar::String = "standard"
    starttime::Union{Nothing, DateTime} = nothing # default from forcing netCDF
    endtime::Union{Nothing, DateTime} = nothing   # default from forcing netCDF
    const time_units::String = CFTime.DEFAULT_TIME_UNITS
    timestepsecs::Union{Nothing, Float64} = nothing # default from forcing netCDF
end

# Logging related configurations
@option struct LoggingSection <: AbstractConfigSection
    silent::Bool = false
    loglevel::String = "info"
    path_log::String = "log.txt"
end

# Water demand configurations
@option struct WaterDemandSubSection <: AbstractConfigSection
    domestic__flag::Bool = false
    industry__flag::Bool = false
    livestock__flag::Bool = false
    paddy__flag::Bool = false
    nonpaddy__flag::Bool = false
end

# Model configurations
@option struct ModelSection <: AbstractConfigSection
    # General
    type::String
    cold_start__flag::Bool = true
    cell_length_in_meter__flag::Bool = false
    reservoir__flag::Bool = false
    # Model types sbm and sbm_gwf
    snow__flag::Bool = false
    snow_gravitational_transport__flag::Bool = false
    glacier__flag::Bool = false
    soil_infiltration_reduction__flag::Bool = false
    soil_layer__thickness::Vector{Int} = [100, 300, 800]
    saturated_hydraulic_conductivity_profile::String = "exponential"
    # Routing method
    land_routing::String = "kinematic-wave"
    river_routing::String = "kinematic-wave"
    # Kinematic wave routing
    pit__flag::Bool = false
    river_streamorder__min_count::Int = 6
    land_streamorder__min_count::Int = 5
    kinematic_wave__adaptive_time_step_flag::Bool = false
    river_kinematic_wave__time_step::Float64 = 900.0
    land_kinematic_wave__time_step::Float64 = 3600.0
    # Local inertial routing
    river_local_inertial_flow__alpha_coefficient::Float64 = 0.7
    land_local_inertial_flow__alpha_coefficient::Float64 = 0.7
    land_local_inertial_flow__theta_coefficient::Float64 = 0.8
    river_water_flow_threshold__depth = 1e-3
    land_surface_water_flow_threshold__depth = 1e-3
    river_water_flow__froude_limit_flag = true
    land_surface_water_flow__froude_limit_flag = true
    floodplain_1d__flag::Bool = false
    # Groundwater flow
    conductivity_profile::String = "uniform"
    drain__flag::Bool = false
    constanthead__flag::Bool = false
    subsurface_water_flow__alpha_coefficient::Float64 = 0.25
    # Model type sediment
    rainfall_erosion::String = "answers"
    overland_flow_erosion::String = "answers"
    run_river_model__flag::Bool = false
    land_transport::String = "yalin"
    river_transport::String = "bagnold"
    # Water demand
    water_demand::WaterDemandSubSection = WaterDemandSubSection()
end

# State configurations
@option struct StateSection <: AbstractConfigSection
    path_input::Union{Nothing, String} = nothing
    path_output::Union{Nothing, String} = nothing
    # Variable name mapping
    variables::MaybePropertyDict = Dict{String, Any}()
end

# Input configurations
@option struct InputSection <: AbstractConfigSection
    # Flow direction and modelling domains
    path_forcing::String
    path_static::String
    basin__local_drain_direction::String
    basin_pit_location__mask::Union{Nothing, String} = nothing
    river_location__mask::String
    reservoir_area__count::String = ""
    reservoir_location__count::String = ""
    subbasin_location__count::String
    # Variable name mappings
    forcing::MaybePropertyDict
    static::MaybePropertyDict
    cyclic::MaybePropertyDict = Dict{String, Any}()
end

# API variable configuration
@option struct APISection <: AbstractConfigSection
    variables::Vector{String} = []
end

# Whether certain relevant sections exist in the TOML
@option struct HasSection <: AbstractConfigSection
    API::Bool
    model_water_demand::Bool
    input_cyclic::Bool
end

# Fields with a default value are optional
@option struct Config
    dir_input::Union{Nothing, String} = nothing
    dir_output::Union{Nothing, String} = nothing
    fews_run__flag::Bool = false
    time::TimeSection = TimeSection()
    logging::LoggingSection = LoggingSection()
    model::ModelSection
    state::StateSection = StateSection()
    input::InputSection
    output::MaybePropertyDict = Dict{String, Any}()
    API::APISection = APISection()
    has_section::HasSection
    path::Union{Nothing, String}  # path to the TOML file, or nothing
end

# Pretty printing the config
function Base.show(io::IO, c::Config)
    println(io, "Wflow Config")
    for field in fieldnames(typeof(c))
        f = getfield(c, field)
        f === nothing || println(io, "\t$field\t= $f")
    end
end

# Create a nested PropertyDict to retrieve data from a nested dict as
# if it were a nested struct
function nested_property_dict(dict::AbstractDict)
    new_dict = Dict{String, Any}()
    for (key, value) in dict
        if value isa AbstractDict
            new_dict[key] = nested_property_dict(value)
        else
            new_dict[key] = value
        end
    end

    return PropertyDict(new_dict)
end

"""
    Config(path::AbstractString; override...)
    Config(dict::AbstractDict; path::Union{Nothing, String} = nothing, override...)

Struct that contains the parsed TOML configurations, as well as a reference to the TOML path
if it exists. The object behaves largely like an immutable nested struct, apart from the flexible
sections (see the fields with `MaybePropertyDict` type in `config.jl`) which are nested property
dicts (see PropertyDicts.jl) and thus mutable.
Values from the TOML (or dict) kan be overridden by passing keyword arguments of the form
`path_to_value_with_underscores = value`, for example `Config(path; logging_silent = true)`.
"""
function Config(path::AbstractString; override...)
    dict = TOML.parsefile(path)
    Config(dict; path, override...)
end

function Config(dict::AbstractDict; path::Union{Nothing, String} = nothing, override...)
    # These sections are particularly checked here because they are mandatory and
    # required to initialize the has_section section
    @assert haskey(dict, "model") "The TOML is missing the mandatory [model] section."
    @assert haskey(dict, "input") "The TOML is missing the mandatory [input] section."
    # When optional values are omitted they obtain the default value, and it cannot be derived
    # whether this value was obtained from the TOML or not. Therefore this information is gathered
    # here.
    dict["has_section"] = Dict(
        "API" => haskey(dict, "API"),
        "model_water_demand" => haskey(dict["model"], "water_demand"),
        "input_cyclic" => haskey(dict["input"], "cyclic"),
    )
    dict["path"] = path
    config = from_dict(Config, dict; override...)

    # Flexible input
    input = config.input
    @reset input.forcing = nested_property_dict(input.forcing)
    @reset input.static = nested_property_dict(input.static)
    @reset input.cyclic = nested_property_dict(input.cyclic)
    @reset config.input = input

    # Flexible output
    @reset config.output = nested_property_dict(config.output)

    # Flexible state
    state = config.state
    @reset state.variables = nested_property_dict(state.variables)
    @reset config.state = state
end

Base.dirname(config::Config) = isnothing(config.path) ? nothing : dirname(config.path)

"Construct a path relative to both the TOML directory and the optional `dir_input`"
input_path(config::Config, path::Union{Nothing, AbstractString}) =
    normpath(dirname(config), config.dir_input, path)

"Construct a path relative to both the TOML directory and the optional `dir_output`"
output_path(config::Config, path::Union{Nothing, AbstractString}) =
    normpath(dirname(config), config.dir_output, path)