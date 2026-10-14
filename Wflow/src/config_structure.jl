#=
Code to support input and output of data and configuration.
Input data can be loaded from netCDF files.
Output data can be written to netCDF or CSV files.
For configuration files we use TOML.
=#

# Option enumerators
@enumx RoutingType kinematic_wave manning_staggered local_inertial
@enumx ModelType sbm sbm_gwf sediment
@enumx CalendarType standard gregorian proleptic_gregorian julian noleap _365_day all_leap _366_day _360_day
@enumx GwfConductivityProfileType uniform exponential
@enumx VerticalConductivityProfile exponential exponential_constant layered layered_exponential
@enumx RainfallErosionType answers eurosem
@enumx OverlandFlowErosionType answers
@enumx LandTransportType yalinpart govers yalin
@enumx RiverTransportType bagnold engelund yang kodatie molinas
@enumx ReducerType maximum minimum mean median sum first last only

abstract type AbstractConfigSection end
const PropertyDictType = PropertyDict{String, Any, Dict{String, Any}}

"""
    @config_section mutable struct Section <: AbstractConfigSection ... end

Like `Base.@kwdef`, and additionally define a type stable
`init_config_section_default(::Type{Section}, dict::AbstractDict{String})` that parses a
TOML dictionary field by field, using the declared field types and defaults. This avoids
splatting a `Dict{Symbol, Any}` into the keyword constructor, which needs dynamic dispatch
and is therefore incompatible with `juliac --trim`.
"""
macro config_section(expr)
    @assert Meta.isexpr(expr, :struct)
    type_name = expr.args[2]
    section_type = Meta.isexpr(type_name, :<:) ? type_name.args[1] : type_name
    fields = Tuple{Symbol, Any, Any}[]
    for arg in expr.args[3].args
        arg isa LineNumberNode && continue
        has_default = Meta.isexpr(arg, :(=))
        lhs = has_default ? arg.args[1] : arg
        field_name, field_type =
            Meta.isexpr(lhs, :(::)) ? (lhs.args[1], lhs.args[2]) : (lhs, :Any)
        key = String(field_name)
        fallback = has_default ? arg.args[2] : :(missing_field_error($section_type, $key))
        push!(fields, (field_name, field_type, fallback))
    end
    assignments = map(fields) do (field_name, field_type, fallback)
        key = String(field_name)
        :(
            $field_name =
                haskey(dict, $key) ?
                convert_field($section_type, $key, $field_type, dict[$key]) : $fallback
        )
    end
    field_names = first.(fields)
    return esc(
        quote
            Base.@kwdef $expr
            function init_config_section_default(
                    ::Type{$section_type},
                    dict::AbstractDict{String},
                )::$section_type
                warn_unknown_keys($section_type, dict)
                $(assignments...)
                return $section_type($(field_names...))
            end
        end,
    )
end

###
### Time section
###

# Time related configurations
@config_section mutable struct TimeSection <: AbstractConfigSection
    calendar::CalendarType.T = CalendarType.standard
    starttime::Union{Nothing, DateTime} = nothing # default from forcing netCDF
    endtime::Union{Nothing, DateTime} = nothing   # default from forcing netCDF
    time_units::String = CFTime.DEFAULT_TIME_UNITS
    timestepsecs::Union{Nothing, Float64} = nothing # default from forcing netCDF
end

###
### Logging section
###

# Logging related configurations
@config_section mutable struct LoggingSection <: AbstractConfigSection
    silent::Bool = false
    loglevel::LogLevel = Info
    path_log::String = "log.txt"
end

###
### Model section
###

# Water demand configurations
@config_section mutable struct WaterDemandSubSection <: AbstractConfigSection
    domestic__flag::Bool = false
    industry__flag::Bool = false
    livestock__flag::Bool = false
    paddy__flag::Bool = false
    nonpaddy__flag::Bool = false
    _was_specified::Bool = true
end

# Model configurations
@config_section mutable struct ModelSection <: AbstractConfigSection
    # General
    type::ModelType.T
    cold_start__flag::Bool = true
    cell_length_in_meter__flag::Bool = false
    reservoir__flag::Bool = false
    # Model types sbm and sbm_gwf
    snow__flag::Bool = false
    snow_gravitational_transport__flag::Bool = false
    glacier__flag::Bool = false
    soil_infiltration_reduction__flag::Bool = false
    soil_layer__thickness::Vector{Int} = [100, 300, 800]
    saturated_hydraulic_conductivity_profile::VerticalConductivityProfile.T =
        VerticalConductivityProfile.exponential
    water_mass_balance__flag::Bool = false
    # Land surface temperature
    land_surface_temperature__flag::Bool = false
    land_surface_wind__speed_reference_height::Float64 = 100.0
    # Routing method
    land_routing::RoutingType.T = RoutingType.kinematic_wave
    river_routing::RoutingType.T = RoutingType.kinematic_wave
    # Kinematic wave routing
    pit__flag::Bool = false
    river_streamorder__min_count::Int = 6
    land_streamorder__min_count::Int = 5
    kinematic_wave__adaptive_time_step_flag::Bool = false
    river_kinematic_wave__time_step::Float64 = 900.0
    land_kinematic_wave__time_step::Float64 = 3600.0
    subsurface_kinematic_wave__time_step::Float64 = 86400.0
    subsurface_kinematic_wave__alpha_coefficient::Float64 = 1.0
    river_subsurface_exchange_head_based__flag::Bool = false
    # Local inertial routing
    river_local_inertial_flow__alpha_coefficient::Float64 = 0.7
    land_local_inertial_flow__alpha_coefficient::Float64 = 0.7
    land_local_inertial_flow__theta_coefficient::Float64 = 1.0
    land_surface_water_flow_threshold__depth::Float64 = 1.0e-3
    river_water_flow__froude_limit_flag::Bool = true
    land_surface_water_flow__froude_limit_flag::Bool = true
    # Local inertial and staggered manning river flow routing
    river_water_flow_threshold__depth::Float64 = 1.0e-3
    floodplain_1d__flag::Bool = false
    # Staggered manning river flow routing
    river_staggered_manning_flow__alpha_coefficient::Float64 = 0.7
    # Groundwater flow
    conductivity_profile::GwfConductivityProfileType.T = GwfConductivityProfileType.uniform
    drain__flag::Bool = false
    constanthead__flag::Bool = false
    subsurface_water_flow__alpha_coefficient::Float64 = 0.25
    # Model type sediment/
    rainfall_erosion::RainfallErosionType.T = RainfallErosionType.answers
    overland_flow_erosion::OverlandFlowErosionType.T = OverlandFlowErosionType.answers
    run_river_model__flag::Bool = false
    land_transport::LandTransportType.T = LandTransportType.yalin
    river_transport::RiverTransportType.T = RiverTransportType.bagnold
    # Water demand
    water_demand::WaterDemandSubSection = WaterDemandSubSection(; _was_specified = false)
end

###
### State section
###

# State configurations
@config_section mutable struct StateSection <: AbstractConfigSection
    path_input::Union{Nothing, String} = nothing
    path_output::Union{Nothing, String} = nothing
    # Variable name mapping
    variables::PropertyDictType = PropertyDict(Dict{String, Any}())
end

###
### Input section
###

@config_section mutable struct InputEntry <: AbstractConfigSection
    # Option 1
    netcdf_variable_name::Union{Nothing, String} = nothing
    scale::Vector{Float64} = [1.0]
    _do_scaling::Bool = !all(isone, scale)
    _scale_scalar::Bool = isone(length(scale))
    offset::Vector{Float64} = [0.0]
    _do_offsetting::Bool = !all(iszero, offset)
    _offset_scalar::Bool = isone(length(offset))
    layer::Union{Nothing, Vector{Int}} = nothing
    # Option 2
    value::Union{Nothing, Vector{Float64}} = nothing
    # Option 3
    external_name::Union{Nothing, String} = nothing
end

variable_name(var::InputEntry) =
    isnothing(var.netcdf_variable_name) ? var.external_name : var.netcdf_variable_name

@config_section mutable struct InputEntries <: AbstractConfigSection
    dict::Dict{String, InputEntry} = Dict{String, InputEntry}()
end

Base.haskey(input_entries::InputEntries, key) = haskey(input_entries.dict, key)
Base.getindex(input_entries::InputEntries, key) = input_entries.dict[key]
Base.setindex!(input_entries::InputEntries, value, key) = (input_entries.dict[key] = value)
Base.setindex!(input_entries::InputEntries, value::AbstractDict, key) =
    (input_entries.dict[key] = init_config_section(InputEntry, value))
Base.keys(input_entries::InputEntries) = keys(input_entries.dict)
Base.iterate(input_entries::InputEntries) = iterate(input_entries.dict)
Base.iterate(input_entries::InputEntries, state) = iterate(input_entries.dict, state)

# Input configurations
@config_section mutable struct InputSection <: AbstractConfigSection
    # Flow direction and modelling domains
    path_forcing::String
    path_static::String
    basin__local_drain_direction::String
    basin_pit_location__mask::Union{Nothing, String} = nothing
    river_location__mask::String
    reservoir_area__count::Union{Nothing, String} = nothing
    reservoir_location__count::Union{Nothing, String} = nothing
    subbasin_location__count::String
    subbasin_active_location__count::Union{Nothing, Vector{Int}} = nothing
    # Variable name mappings
    forcing::InputEntries
    static::InputEntries
    cyclic::InputEntries = InputEntries()
    _location_maps::PropertyDictType
end

const input_field_names = String.(fieldnames(InputSection))

###
### Output section
###

@config_section mutable struct CoordinateSection <: AbstractConfigSection
    x::Float64
    y::Float64
end

@config_section mutable struct IndexSection <: AbstractConfigSection
    i::Union{Nothing, Int} = nothing
    x::Union{Nothing, Int} = nothing
    y::Union{Nothing, Int} = nothing
    _was_specified::Bool = true
end

@config_section mutable struct CSVColumn <: AbstractConfigSection
    header::String
    parameter::String
    reducer::Union{Nothing, ReducerType.T} = nothing
    layer::Union{Nothing, Int} = nothing
    index::IndexSection = IndexSection(; _was_specified = false)
    map::Union{Nothing, String} = nothing
    coordinate::Union{Nothing, CoordinateSection} = nothing
end

@config_section mutable struct CSVSection <: AbstractConfigSection
    path::String
    column::Vector{CSVColumn} = CSVColumn[]
    _was_specified::Bool = true
end

@config_section mutable struct NetCDFScalarVariable <: AbstractConfigSection
    name::String
    parameter::String
    reducer::Union{Nothing, ReducerType.T} = nothing
    layer::Union{Nothing, Int} = nothing
    index::IndexSection = IndexSection(; _was_specified = false)
    map::Union{Nothing, String} = nothing
    coordinate::Union{Nothing, CoordinateSection} = nothing
    location::Union{Nothing, String} = nothing
    _location_dim::String
end

@config_section mutable struct NetCDFScalarSection <: AbstractConfigSection
    path::String
    variable::Vector{NetCDFScalarVariable} = NetCDFScalarVariable[]
    _was_specified::Bool = true
end

@config_section mutable struct NetCDFGridSection <: AbstractConfigSection
    path::String
    compressionlevel::Int = 0
    variables::PropertyDictType = PropertyDict(Dict{String, Any}())
    _was_specified::Bool = true
end

@config_section mutable struct OutputSection <: AbstractConfigSection
    csv::CSVSection = CSVSection(; _was_specified = false, path = "")
    netcdf_scalar::NetCDFScalarSection =
        NetCDFScalarSection(; _was_specified = false, path = "")
    netcdf_grid::NetCDFGridSection = NetCDFGridSection(; _was_specified = false, path = "")
end

###
### API section
###

# API variable configuration
@config_section mutable struct APISection <: AbstractConfigSection
    variables::Vector{String} = String[]
    _was_specified::Bool = true
end

###
### Config section (top level)
###

# Fields with a default value are optional
@config_section mutable struct Config <: AbstractConfigSection
    wflow_version::Union{Nothing, VersionNumber} = nothing
    dir_input::String = "."
    dir_output::String = "."
    fews_run__flag::Bool = false
    time::TimeSection = TimeSection()
    logging::LoggingSection = LoggingSection()
    model::ModelSection
    state::StateSection = StateSection()
    input::InputSection
    output::OutputSection = OutputSection()
    API::APISection = APISection(; _was_specified = false)
    path::String  # path to the TOML file
end
