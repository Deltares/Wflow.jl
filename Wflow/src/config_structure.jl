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

abstract type AbstractConfigSection end
const PropertyDictType = PropertyDict{String, Any, Dict{String, Any}}

###
### Time section
###

# Time related configurations
@kwdef mutable struct TimeSection <: AbstractConfigSection
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
@kwdef struct LoggingSection <: AbstractConfigSection
    silent::Bool = false
    loglevel::LogLevel = Info
    path_log::String = "log.txt"
end

###
### Model section
###

# Water demand configurations
@kwdef struct WaterDemandSubSection <: AbstractConfigSection
    domestic__flag::Bool = false
    industry__flag::Bool = false
    livestock__flag::Bool = false
    paddy__flag::Bool = false
    nonpaddy__flag::Bool = false
    _was_specified::Bool = true
end

# Model configurations
@kwdef mutable struct ModelSection <: AbstractConfigSection
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
@kwdef struct StateSection <: AbstractConfigSection
    path_input::Union{Nothing, String} = nothing
    path_output::Union{Nothing, String} = nothing
    # Variable name mapping
    variables::PropertyDictType = PropertyDict(Dict{String, Any}())
end

###
### Input section
###

@kwdef struct InputEntry <: AbstractConfigSection
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

@kwdef struct InputEntries <: AbstractConfigSection
    dict::Dict{String, InputEntry} = Dict()
end

Base.haskey(input_entries::InputEntries, key) = haskey(input_entries.dict, key)
Base.getindex(input_entries::InputEntries, key) = input_entries.dict[key]
Base.keys(input_entries::InputEntries) = keys(input_entries.dict)
Base.iterate(input_entries::InputEntries) = iterate(input_entries.dict)
Base.iterate(input_entries::InputEntries, state) = iterate(input_entries.dict, state)

# Input configurations
@kwdef struct InputSection <: AbstractConfigSection
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

###
### Output section
###

@kwdef struct CoordinateSection <: AbstractConfigSection
    x::Float64
    y::Float64
end

@kwdef struct IndexSection <: AbstractConfigSection
    i::Union{Nothing, Int} = nothing
    x::Union{Nothing, Int} = nothing
    y::Union{Nothing, Int} = nothing
    _was_specified::Bool = true
end

@kwdef struct CSVColumn <: AbstractConfigSection
    header::String
    parameter::String
    reducer::Union{Nothing, ReducerType.T} = nothing
    layer::Union{Nothing, Int} = nothing
    index::IndexSection = IndexSection(; _was_specified = false)
    map::Union{Nothing, String} = nothing
    coordinate::Union{Nothing, CoordinateSection} = nothing
end

@kwdef struct CSVSection <: AbstractConfigSection
    path::String
    column::Vector{CSVColumn} = []
    _was_specified::Bool = true
end

@kwdef struct NetCDFScalarVariable <: AbstractConfigSection
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

@kwdef struct NetCDFScalarSection <: AbstractConfigSection
    path::String
    variable::Vector{NetCDFScalarVariable} = []
    _was_specified::Bool = true
end

@kwdef struct NetCDFGridSection <: AbstractConfigSection
    path::String
    compressionlevel::Int = 0
    variables::PropertyDictType = PropertyDict(Dict{String, Any}())
    _was_specified::Bool = true
end

@kwdef struct OutputSection <: AbstractConfigSection
    csv::CSVSection = CSVSection(; _was_specified = false, path = "")
    netcdf_scalar::NetCDFScalarSection =
        NetCDFScalarSection(; _was_specified = false, path = "")
    netcdf_grid::NetCDFGridSection = NetCDFGridSection(; _was_specified = false, path = "")
end

###
### API section
###

# API variable configuration
@kwdef struct APISection <: AbstractConfigSection
    variables::Vector{String} = []
    _was_specified::Bool = true
end

###
### Config section (top level)
###

# Fields with a default value are optional
@kwdef struct Config <: AbstractConfigSection
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