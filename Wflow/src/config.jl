#=
Code to support input and output of data and configuration.
Input data can be loaded from netCDF files.
Output data can be written to netCDF or CSV files.
For configuration files we use TOML.
=#

# Option enumerators
@enumx RoutingType kinematic_wave local_inertial
@enumx ModelType sbm sbm_gwf sediment
@enumx CalendarType standard gregorian proleptic_gregorian julian noleap _365_day all_leap _366_day _360_day
@enumx ConductivityProfileType uniform exponential
@enumx SHCPType exponential exponential_constant layered layered_exponential
@enumx RainfallErosionType answers eurosem
@enumx OverlandFlowErosionType answers
@enumx LandTransportType yalinpart govers yalin
@enumx RiverTransportType bagnold engelund yang kodatie molinas
@enumx ReducerType maximum minimum mean median sum first last only

Base.convert(::Type{String}, calendar_type::CalendarType.T) =
    remove_leading_underscore(String(Symbol(calendar_type)))

abstract type AbstractConfigSection end

# e.g. section_name(InputSection) = "[input]"
function section_name(::Type{<:AbstractConfigSection})
    name = lowercase(replace(String(nameof(InputSection)), "Section" => ""))
    return "[$name]"
end

# Don't error on extra fields that aren't hardcoded below
Configurations.ignore_extra(::Type{<:AbstractConfigSection}) = true

const PropertyDictType = PropertyDict{String, Any, Dict{String, Any}}

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

# Let Configurations.jl wrap dict fields in a PropertyDict
Configurations.from_dict(
    ::Type{<:AbstractConfigSection},
    ::Type{PropertyDictType},
    dict::AbstractDict{String},
) = nested_property_dict(dict)

# Handle options starting with a number
add_leading_underscore(s::String) = isdigit(first(s)) ? "_$s" : s
remove_leading_underscore(s::String) = (first(s) == '_') ? s[2:end] : s

# Get printable tuple of possible options
option_names(::Type{T}) where {T <: EnumX.Enum} =
    remove_leading_underscore.(String.(Symbol.(instances(T))))

# Account for the dashes in the routing type names, which cannot be present in the enum
option_names(::Type{RoutingType.T}) = ("kinematic-wave", "local-inertial")

# Let Configurations.jl automatically convert strings into enumerator instances
function Configurations.from_dict(
    S::Type{<:AbstractConfigSection},
    ::OptionField{NAME},
    ::Type{T},
    option_string::String,
) where {T <: Union{Nothing, EnumX.Enum}, NAME}
    option = Symbol(add_leading_underscore(replace(option_string, "-" => "_")))
    for instance in instances(T)
        if option == Symbol(instance)
            return instance
        end
    end
    throw(
        ArgumentError(
            "Invalid value $NAME = $option_string in the $(section_name(S)) section of the TOML, must be one of $(option_names(T)).",
        ),
    )
end

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

###
### Time section
###

# Time related configurations
@option mutable struct TimeSection <: AbstractConfigSection
    const calendar::CalendarType.T = CalendarType.standard
    starttime::Union{Nothing, DateTime} = nothing # default from forcing netCDF
    endtime::Union{Nothing, DateTime} = nothing   # default from forcing netCDF
    const time_units::String = CFTime.DEFAULT_TIME_UNITS
    timestepsecs::Union{Nothing, Float64} = nothing # default from forcing netCDF
end

###
### Logging section
###

# Logging related configurations
@option struct LoggingSection <: AbstractConfigSection
    silent::Bool = false
    loglevel::LogLevel = Info
    path_log::String = "log.txt"
end

const log_level_map::Dict{Union{Int, String}, LogLevel} =
    Dict("debug" => Debug, "info" => Info, "warn" => Warn, "error" => Error)

parse_loglevel(level::Int) = LogLevel(level)
parse_loglevel(level::String) = get(log_level_map, lowercase(level), missing)

function Configurations.from_dict(
    ::Type{LoggingSection},
    ::Type{LogLevel},
    input_string::String,
)
    level = parse_loglevel(input_string)
    ismissing(level) && throw(
        ArgumentError(
            "Invalid loglevel = $input_string in the $(section_name(LoggingSection)) section of the TOML, must be one of $(keys(log_level_map)).",
        ),
    )
    return level
end

###
### Model section
###

# Water demand configurations
@option struct WaterDemandSubSection <: AbstractConfigSection
    domestic__flag::Bool = false
    industry__flag::Bool = false
    livestock__flag::Bool = false
    paddy__flag::Bool = false
    nonpaddy__flag::Bool = false
    _was_specified::Bool = true
end

# Model configurations
@option mutable struct ModelSection <: AbstractConfigSection
    # General
    const type::ModelType.T
    cold_start__flag::Bool = true
    const cell_length_in_meter__flag::Bool = false
    const reservoir__flag::Bool = false
    # Model types sbm and sbm_gwf
    const snow__flag::Bool = false
    const snow_gravitational_transport__flag::Bool = false
    const glacier__flag::Bool = false
    const soil_infiltration_reduction__flag::Bool = false
    const soil_layer__thickness::Vector{Int} = [100, 300, 800]
    const saturated_hydraulic_conductivity_profile::SHCPType.T = SHCPType.exponential
    # Routing method
    const land_routing::RoutingType.T = RoutingType.kinematic_wave
    const river_routing::RoutingType.T = RoutingType.kinematic_wave
    # Kinematic wave routing
    const pit__flag::Bool = false
    const river_streamorder__min_count::Int = 6
    const land_streamorder__min_count::Int = 5
    const kinematic_wave__adaptive_time_step_flag::Bool = false
    const river_kinematic_wave__time_step::Float64 = 900.0
    const land_kinematic_wave__time_step::Float64 = 3600.0
    # Local inertial routing
    const river_local_inertial_flow__alpha_coefficient::Float64 = 0.7
    const land_local_inertial_flow__alpha_coefficient::Float64 = 0.7
    const land_local_inertial_flow__theta_coefficient::Float64 = 0.8
    const river_water_flow_threshold__depth = 1e-3
    const land_surface_water_flow_threshold__depth = 1e-3
    const river_water_flow__froude_limit_flag = true
    const land_surface_water_flow__froude_limit_flag = true
    const floodplain_1d__flag::Bool = false
    # Groundwater flow
    const conductivity_profile::ConductivityProfileType.T = ConductivityProfileType.uniform
    const drain__flag::Bool = false
    const constanthead__flag::Bool = false
    subsurface_water_flow__alpha_coefficient::Float64 = 0.25
    # Model type sediment/
    const rainfall_erosion::RainfallErosionType.T = RainfallErosionType.answers
    const overland_flow_erosion::OverlandFlowErosionType.T = OverlandFlowErosionType.answers
    const run_river_model__flag::Bool = false
    const land_transport::LandTransportType.T = LandTransportType.yalin
    const river_transport::RiverTransportType.T = RiverTransportType.bagnold
    # Water demand
    const water_demand::WaterDemandSubSection =
        WaterDemandSubSection(; _was_specified = false)
end

###
### State section
###

# State configurations
@option struct StateSection <: AbstractConfigSection
    path_input::Union{Nothing, String} = nothing
    path_output::Union{Nothing, String} = nothing
    # Variable name mapping
    variables::PropertyDictType = PropertyDict(Dict{String, Any}())
end

###
### Input section
###

@option struct InputEntry <: AbstractConfigSection
    # Option 1
    netcdf_variable_name::Union{Nothing, String} = nothing # Comes from dict["netcdf"]["variable"]["name"]
    scale::Vector{Float64} = [1.0]
    offset::Vector{Float64} = [0.0]
    layer::Union{Nothing, Vector{Int}} = nothing
    # Option 2
    value::Any = nothing
    # Option 3
    standard_name::Union{Nothing, String} = nothing
end

variable_name(var::InputEntry) =
    isnothing(var.netcdf_variable_name) ? var.standard_name : var.netcdf_variable_name

# Get a value from a nested dict if each intermediate
# key exists, otherwise return nothing
function get_nested(dict::AbstractDict, keys)
    return if length(keys) == 1
        get(dict, only(keys), nothing)
    else
        key = first(keys)
        if haskey(dict, key) && dict[key] isa AbstractDict
            get_nested(dict[key], view(keys, 2:length(keys)))
        else
            nothing
        end
    end
end

# Option 1 and 2
function Configurations.from_dict(::Type{InputEntry}, dict::AbstractDict{String})
    netcdf_variable_name = get_nested(dict, ["netcdf", "variable", "name"])
    value = get(dict, "value", nothing)

    if !isnothing(netcdf_variable_name)
        dict["netcdf_variable_name"] = netcdf_variable_name
        # Option 1
        len = missing

        for field_name in ("scale", "offset", "layer")
            if haskey(dict, field_name)
                value = dict[field_name]
                # Convert scalar to vector
                value = (value isa Number) ? [value] : value
                dict[field_name] = value
                if !ismissing(len) && length(value) != len
                    throw(
                        ArgumentError(
                            "In the TOML [input] section with netcdf variable name $netcdf_variable_name," *
                            "the 'scale' 'offset' and 'layer' values must all have the same length if provided.",
                        ),
                    )
                else
                    len = length(value)
                end
            end
        end

        # Add a proper amount of the default values
        len = coalesce(len, 0)
        !haskey(dict, "scale") && (dict["scale"] = ones(len))
        !haskey(dict, "offset") && (dict["offset"] = zeros(len))
        has_layer = haskey(dict, "layer")

        for (i, (scale, offset)) in enumerate(zip(dict["scale"], dict["offset"]))
            # Modifier info
            msg = "NetCDF parameter $netcdf_variable_name is modified with scale $scale and offset $offset"
            has_layer && (msg *= "at index $(dict["layer"][i])")
            @info "$msg."
        end

        # Invoke default method
        return invoke(from_dict, Tuple{Type, AbstractDict{String}}, InputEntry, dict)
    elseif !isnothing(value)
        # Option 2
        return InputEntry(; value)
    end
    throw(
        ArgumentError(
            "Couldn't parse the section under [input] in the TOML with data $dict.",
        ),
    )
end

# Option 3 where only a standard name is given
Configurations.from_dict(::Type{InputEntry}, standard_name::String) =
    InputEntry(; standard_name)
Base.convert(::Type{InputEntry}, standard_name::String) = InputEntry(; standard_name)

@option struct InputEntries <: AbstractConfigSection
    dict::Dict{String, InputEntry} = Dict()
end

Base.haskey(input_entries::InputEntries, key) = haskey(input_entries.dict, key)
Base.getindex(input_entries::InputEntries, key) = input_entries.dict[key]
Base.keys(input_entries::InputEntries) = keys(input_entries.dict)
Base.iterate(input_entries::InputEntries) = iterate(input_entries.dict)
Base.iterate(input_entries::InputEntries, state) = iterate(input_entries.dict, state)

function Configurations.from_dict(::Type{InputEntries}, dict::AbstractDict{String})
    InputEntries(Dict(key => from_dict(InputEntry, value) for (key, value) in dict))
end

# Input configurations
@option struct InputSection <: AbstractConfigSection
    # Flow direction and modelling domains
    path_forcing::String
    path_static::String
    basin__local_drain_direction::String
    basin_pit_location__mask::Union{Nothing, String} = nothing
    river_location__mask::String
    reservoir_area__count::Union{Nothing, String} = nothing
    reservoir_location__count::Union{Nothing, String} = nothing
    subbasin_location__count::String
    # Variable name mappings
    forcing::InputEntries
    static::InputEntries
    cyclic::InputEntries = InputEntries()
    flexible::PropertyDictType
end

const input_field_names = String.(fieldnames(Wflow.InputSection))

function Configurations.from_dict(::Type{InputSection}, dict::AbstractDict{String})
    # Move flexible part of the input section into input.flexible
    dict["flexible"] = filter(kv -> kv[1] ∉ input_field_names, dict)
    # Invoke default method
    return invoke(from_dict, Tuple{Type, AbstractDict{String}}, InputSection, dict)
end

###
### Output section
###

@option struct CoordinateSection <: AbstractConfigSection
    x::Float64
    y::Float64
end

@option struct IndexSection <: AbstractConfigSection
    x::Int
    y::Int
end

# TODO: Validate on construction
@option struct CSVColumn <: AbstractConfigSection
    header::String
    parameter::String
    reducer::Union{Nothing, ReducerType.T} = nothing
    layer::Union{Nothing, Int} = nothing
    index::Union{Nothing, Int, IndexSection} = nothing
    map::Union{Nothing, String} = nothing
    coordinate::Union{Nothing, CoordinateSection} = nothing
end

@option struct CSVSection <: AbstractConfigSection
    path::String
    column::Vector{CSVColumn} = []
    _was_specified::Bool = true
end

@option struct NetCDFScalarVariable <: AbstractConfigSection
    name::String
    parameter::String
    reducer::Union{Nothing, ReducerType.T} = nothing
    layer::Union{Nothing, Int} = nothing
    index::Union{Nothing, Int, IndexSection} = nothing
    map::Union{Nothing, String} = nothing
    coordinate::Union{Nothing, CoordinateSection} = nothing
    location::Union{Nothing, String} = nothing
    _location_dim::String
end

@option struct NetCDFScalarSection <: AbstractConfigSection
    path::String
    variable::Vector{NetCDFScalarVariable} = []
    _was_specified::Bool = true
end

function Configurations.from_dict(::Type{NetCDFScalarVariable}, dict::AbstractDict{String})
    for key in ("name", "parameter")
        !haskey(dict, key) && throw(
            ArgumentError(
                "A variable in a [[output.netcdf_scalar.variable]] section in the TOML is missing the mandatory field $key.",
            ),
        )
    end

    if haskey(dict, "map")
        dict["_location_dim"] = string(dict["name"], '_', dict["map"])
    elseif haskey(dict, "location")
        dict["_location_dim"] = dict["location"]
    else
        throw(
            ArgumentError(
                "The variable in the [[output.netcdf_scalar.variable]] section with name = '$name'" *
                " is missing either a 'location' or a 'map' field.",
            ),
        )
    end
    # Invoke default method
    return invoke(from_dict, Tuple{Type, AbstractDict{String}}, NetCDFScalarVariable, dict)
end

@option struct NetCDFGridSection <: AbstractConfigSection
    path::String
    compressionlevel::Int = 0
    variables::PropertyDictType = PropertyDict(Dict{String, Any}())
    _was_specified::Bool = true
end

@option struct OutputSection <: AbstractConfigSection
    csv::CSVSection = CSVSection(; _was_specified = false, path = "")
    netcdf_scalar::NetCDFScalarSection =
        NetCDFScalarSection(; _was_specified = false, path = "")
    netcdf_grid::NetCDFGridSection = NetCDFGridSection(; _was_specified = false, path = "")
end

###
### API section
###

# API variable configuration
@option struct APISection <: AbstractConfigSection
    variables::Vector{String} = []
    _was_specified::Bool = true
end

###
### Config section (top level)
###

# Fields with a default value are optional
@option struct Config <: AbstractConfigSection
    dir_input::Union{Nothing, String} = nothing
    dir_output::Union{Nothing, String} = nothing
    fews_run__flag::Bool = false
    time::TimeSection = TimeSection()
    logging::LoggingSection = LoggingSection()
    model::ModelSection
    state::StateSection = StateSection()
    input::InputSection
    output::OutputSection = OutputSection()
    API::APISection = APISection(; _was_specified = false)
    path::Union{Nothing, String}  # path to the TOML file, or nothing
end

do_water_demand(config::Config) = config.model.water_demand._was_specified
do_api(config::Config) = config.API._was_specified
do_csv(config::Config) = config.output.csv._was_specified
do_netcdf_scalar(config::Config) = config.output.netcdf_scalar._was_specified
do_netcdf_grid(config::Config) = config.output.netcdf_grid._was_specified
do_cyclic(config::Config) = !haskey(config.input.cyclic, "_was_not_specified")

subsurface_routing(config::Config) =
    config.model.type == ModelType.sbm ? RoutingType.kinematic_wave : nothing

# Pretty printing the config
function Base.show(io::IO, c::Config)
    println(io, "Wflow Config")
    for field in fieldnames(typeof(c))
        f = getfield(c, field)
        f === nothing || println(io, "\t$field\t= $f")
    end
end

"""
    Config(path::AbstractString; validate::Bool = true, override...)
    Config(dict::AbstractDict; path::Union{Nothing, String} = nothing, validate::Bool = true, override...)

Struct that contains the parsed TOML configurations, as well as a reference to the TOML path
if it exists. The object behaves largely like an immutable nested struct, apart from the flexible
sections (see the fields with `PropertyDictType` type in `config.jl`) which are nested property
dicts (see PropertyDicts.jl) and thus mutable.
Values from the TOML (or dict) kan be overridden by passing keyword arguments of the form
`path_to_value_with_underscores = value`, for example `Config(path; logging_silent = true)`.
Fields whose name start with an underscore should not be specified in the TOML.
"""
function Config(path::AbstractString; kwargs...)
    dict = TOML.parsefile(path)
    Config(dict; path, kwargs...)
end

function Config(
    dict::AbstractDict;
    path::Union{Nothing, String} = nothing,
    validate::Bool = true,
    override...,
)
    # Add path to config
    dict["path"] = path
    config = from_dict(Config, dict; override...)
    validate && validate_config(config)
    return config
end

Base.dirname(config::Config) = isnothing(config.path) ? nothing : dirname(config.path)

"Construct a path relative to both the TOML directory and the optional `dir_input`"
input_path(config::Config, path::Union{Nothing, AbstractString}) =
    normpath(dirname(config), config.dir_input, path)

"Construct a path relative to both the TOML directory and the optional `dir_output`"
output_path(config::Config, path::Union{Nothing, AbstractString}) =
    normpath(dirname(config), config.dir_output, path)

function validate_config(config::Config)::Nothing
    valid = true

    # check if there is overlap in forcing and cyclic parameters
    if do_cyclic(config)
        overlap = intersect(keys(config.input.forcing), keys(config.input.cyclic))
        if !isempty(overlap)
            @error "These parameters were specified in both the forcing and cyclic" *
                   "sections of the [input] in the TOML: $overlap."
            valid = false
        end
    end

    !valid && error("Invalid TOML file.")
    return nothing
end

# config value converter to write the Config instance to a TOML file
# usage: Configurations.to_toml(custom_convert, path, config)
custom_convert(x::Any) = x
custom_convert(x::EnumX.Enum) = remove_leading_underscore(String(Symbol(x)))
custom_convert(x::RoutingType.T) =
    (x == RoutingType.kinematic_wave) ? "kinematic-wave" : "local-inertial"
custom_convert(x::LogLevel) = lowercase(string(x))
custom_convert(var::InputEntry) = to_dict(var)
custom_convert(::Nothing) = ""