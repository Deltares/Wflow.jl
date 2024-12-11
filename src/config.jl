module config

using Configurations: @option, from_toml
using Dates: DateTime

export Config

abstract type TableOption end

@option struct Model <: TableOption
    type::String = "sbm"
    river_routing::String = "kinematic-wave"
    lakes::Bool = false
    masswasting::Bool = false
    snow::Bool = false
    glacier::Bool = false
    reinit::Bool = true
    reservoirs::Bool = false
    kin_wave_iteration::Bool = false
    thicknesslayers::Vector{Float64} = Float64[]
    min_streamorder_river::Int = 6
    min_streamorder_land::Int = 5
end

@option struct StateVerticalInterceptionVariables <: TableOption
    canopy_storage::String = "canopystorage"
end

@option struct StateVerticalInterception <: TableOption
    variables::StateVerticalInterceptionVariables = StateVerticalInterceptionVariables()
end

@option struct StateVerticalSoilVariables <: TableOption
    satwaterdepth::String = "satwaterdepth"
    tsoil::String = "tsoil"
    ustorelayerdepth::String = "ustorelayerdepth"
end

@option struct StateVerticalSoil <: TableOption
    variables::StateVerticalSoilVariables = StateVerticalSoilVariables()
end

@option struct StateVerticalGlacierVariables <: TableOption
    glacier_store::String = "glacierstore"
end

@option struct StateVerticalGlacier <: TableOption
    variables::StateVerticalGlacierVariables = StateVerticalGlacierVariables()
end

@option struct StateVerticalSnowVariables <: TableOption
    snow_storage::String = "snow"
    snow_water = "snowwater"
end

@option struct StateVerticalSnow <: TableOption
    variables::StateVerticalSnowVariables = StateVerticalSnowVariables()
end

@option struct StateVertical <: TableOption
    satwaterdept::String = "satwaterdepth"
    tsoil::String = "tsoil"
    ustorelayerdepth::String = "ustorelayerdepth"
    snowwater::String = "snowwater"
    canopystorage::String = "canopystorage"
    interception::StateVerticalInterception = StateVerticalInterception()
    soil::StateVerticalSoil = StateVerticalSoil()
    glacier::StateVerticalGlacier = StateVerticalGlacier()
    snow::StateVerticalSnow = StateVerticalSnow()
end

@option struct StateLateralRiverVariables <: TableOption
    q::String = "q_river"
    h::String = "h_river"
    h_av::String = "h_av_river"
end

@option struct StateLateralRiver <: TableOption
    variables::StateLateralRiverVariables = StateLateralRiverVariables()
end

@option struct StateLateralSubsurfaceVariables <: TableOption
    ssf::String = "ssf"
end

@option struct StateLateralSubsurface <: TableOption
    variables::StateLateralSubsurfaceVariables = StateLateralSubsurfaceVariables()
end

@option struct StateLateralLandVariables <: TableOption
    q::String = "q_land"
    h::String = "h_land"
    h_av::String = "h_av_land"
end

@option struct StateLateralLand <: TableOption
    variables::StateLateralLandVariables = StateLateralLandVariables()
end

@option struct StateLateral <: TableOption
    river::StateLateralRiver = StateLateralRiver()
    subsurface::StateLateralSubsurface = StateLateralSubsurface()
    land::StateLateralLand = StateLateralLand()
end

@option struct State <: TableOption
    path_input::String = "."
    path_output::String = "."
    vertical::StateVertical = StateVertical()
    lateral::StateLateral = StateLateral()
end

@option struct InputVerticalInterceptionParameters <: TableOption
    e_r = "EoverR"
end

@option struct InputVerticalInterception <: TableOption
    parameters::InputVerticalInterceptionParameters = InputVerticalInterceptionParameters()
end

@option struct InputVerticalSoilParameters <: TableOption
    alpha_h1::String = "alpha_h1"
    c::String = "c"
    cf_soil::String = "cf_soil"
    f::String = "f"
    h1::String = "h1"
    h2::String = "h2"
    h3_high::String = "h3_high"
    h3_low::String = "h3_low"
    h4::String = "h4"
    infiltcappath::String = "InfiltCapPath"
    infiltcapsoil::String = "InfiltCapSoil"
    theta_r::String = "thetaR"
    theta_s::String = "thetaS"
    maxleakage::String = "MaxLeakage"
    pathfrac::String = "PathFrac"
    rootdistpar::String = "rootdistpar"
    soilthickness::String = "SoilThickness"
    kv_0::String = "KsatVer"
end

@option struct InputVerticalSoil <: TableOption
    parameters::InputVerticalSoilParameters = InputVerticalSoilParameters()
end

@option struct InputVerticalSnowParameters <: TableOption
    tt::String = "TT"
    tti::String = "TTI"
    ttm::String = "TTM"
    cfmax::String = "Cfmax"
end

@option struct InputVerticalSnow <: TableOption
    parameters::InputVerticalSnowParameters = InputVerticalSnowParameters()
end

@option struct InputVerticalGlacierParameters <: TableOption
    glacier_frac::String = "wflow_glacierfrac"
    g_cfmax::String = "G_Cfmax"
    g_ttm::String = "G_TT"
    g_sifrac::String = "G_SIfrac"
end

@option struct InputVerticalGlacierVariables <: TableOption
    glacier_store::String = "wflow_glacierstore"
end

@option struct InputVerticalGlacier <: TableOption
    parameters::InputVerticalGlacierParameters = InputVerticalGlacierParameters()
    variables::InputVerticalGlacierVariables = InputVerticalGlacierVariables()
end

@option struct InputVerticalAtmosphericForcing <: TableOption
    potential_evaporation::String = "pet"
    precipitation::String = "precip"
    temperature::String = "temp"
end

@option struct InputVerticalRunoffParameters <: TableOption
    waterfrac = "waterfrac"
end

@option struct InputVerticalRunoff <: TableOption
    parameters::InputVerticalRunoffParameters = InputVerticalRunoffParameters()
end

@option struct InputVerticalVegetationParameterSet <: TableOption
    leaf_area_index::String = "LAI"
    kext::String = "Kext"
    storage_specific_leaf::String = "Sl"
    storage_wood::String = "Swood"
    rootingdepth::String = "RootingDepth"
    kc::String = "crop_factor"
end

@option struct InputVertical <: TableOption
    altitude::String = "wflow_dem"
    c::String = "c"
    cf_soil::String = "cf_soil"
    cfmax::String = "Cfmax"
    e_r::String = "EoverR"
    infiltcappath::String = "InfiltCapPath"
    infiltcapsoil::String = "InfiltCapSoil"
    kext::String = "Kext"
    kv_0::String = "KsatVer"
    leaf_area_index::String = "LAI"             # Cyclic variable
    m::String = "M"
    maxleakage::String = "MaxLeakage"
    pathfrac::String = "PathFrac"
    potential_evaporation::String = "PET"       # Forcing variable
    precipitation::String = "P"                 # Forcing variable
    rootdistpar::String = "rootdistpar"
    rootingdepth::String = "RootingDepth"
    soilminthickness::String = "SoilMinThickness"
    soilthickness::String = "SoilThickness"
    specific_leaf::String = "Sl"
    storage_wood::String = "Swood"
    temperature::String = "TEMP"                # Forcing variable
    tt::String = "TT"
    tti::String = "TTI"
    ttm::String = "TTM"
    w_soil::String = "wflow_soil"
    water_holding_capacity::String = "WHC"
    waterfrac::String = "WaterFrac"
    theta_r::String = "thetaR"
    theta_s::String = "thetaS"
    interception::InputVerticalInterception = InputVerticalInterception()
    soil::InputVerticalSoil = InputVerticalSoil()
    snow::InputVerticalSnow = InputVerticalSnow()
    glacier::InputVerticalGlacier = InputVerticalGlacier()
    atmospheric_forcing::InputVerticalAtmosphericForcing = InputVerticalAtmosphericForcing()
    runoff::InputVerticalRunoff = InputVerticalRunoff()
    vegetation_parameter_set::InputVerticalVegetationParameterSet =
        InputVerticalVegetationParameterSet()
end

@option struct InputLateralRiver <: TableOption
    length::String = "wflow_riverlength"
    mannings_n::String = "N_River"
    slope::String = "RiverSlope"
    width::String = "wflow_riverwidth"
    bankfull_depth::String = "RiverDepth"
end

@option struct InputLateralSubsurface <: TableOption
    ksathorfrac::String = "KsatHorFrac"
end

@option struct InputLateralLand <: TableOption
    n::String = "N"
    slope::String = "Slope"
    mannings_n::String = "N"
end

@option struct InputLateral <: TableOption
    river::InputLateralRiver = InputLateralRiver()
    subsurface::InputLateralSubsurface = InputLateralSubsurface()
    land::InputLateralLand = InputLateralLand()
end

@option struct Input <: TableOption
    path_forcing::String = "."
    path_static::String = "."
    gauges::String = "wflow_gauges"
    gauges_grdc::String = "wflow_gauges_grdc"
    ldd::String = "wflow_ldd"
    river_location::String = "wflow_river"
    subcatchment::String = "wflow_subcatch"
    forcing::Vector{String} = String[]
    cyclic::Vector{String} = String[]
    vertical::InputVertical = InputVertical()
    lateral::InputLateral = InputLateral()
end

@option struct OutputVertical <: TableOption
    satwaterdepth::String = "satwaterdepth"
    snow::String = "snow"
    tsoil::String = "tsoil"
    ustorelayerdepth::String = "ustorelayerdepth"
    snowwater::String = "snowwater"
    canopystorage::String = "canopystorage"
end

@option struct OutputLateralRiverReservoir <: TableOption
    volume::String = "volume_reservoir"
end

@option struct OutputLateralRiverVariables <: TableOption
    q::String = "q_river"
    q_av::String = "q_river"
    h::String = "h_river"
    reservoir::OutputLateralRiverReservoir = OutputLateralRiverReservoir()
end

@option struct OutputLateralRiver <: TableOption
    variables::OutputLateralRiverVariables = OutputLateralRiverVariables()
end

@option struct OutputLateralSubsurface <: TableOption
    ssf::String = "ssf"
end

@option struct OutputLateralLand <: TableOption
    q::String = "q_land"
    h::String = "h_land"
end

@option struct OutputLateral <: TableOption
    river::OutputLateralRiver = OutputLateralRiver()
    subsurface::OutputLateralSubsurface = OutputLateralSubsurface()
    land::OutputLateralLand = OutputLateralLand()
end

@option struct Output <: TableOption
    path::String = "."
    vertical::OutputVertical = OutputVertical()
    lateral::OutputLateral = OutputLateral()
end

@option struct CSV <: TableOption
    path::String = "."
    column::Vector{Dict{String, String}} = Dict{String, String}[]
end

@option struct Toml <: TableOption
    calendar::String = "standard"
    starttime::DateTime = DateTime(2000, 1, 1)
    endtime::DateTime = DateTime(2000, 2)
    time_units::String = "days since 1900-01-01 00:00:00"
    timestepsecs::Float64 = 86400
    dir_input::String = "data/input"
    dir_output::String = "data/output"
    # Config
    silent::Bool = false
    loglevel::String = "debug"
    path_log::String = "log.txt"
    fews_run::Bool = false
    model::Model = Model()
    state::State = State()
    input::Input = Input()
    output::Output = Output()
    csv::CSV = CSV()
end

struct Config
    toml::Toml
    dir::String
end

Config(toml::Toml) = Config(toml, ".")

"""
    Config(config_path::AbstractString; kwargs...)

Parse a TOML file to a Config. Keys can be overruled using keyword arguments. To overrule
keys from a subsection, e.g. `dt` from the `solver` section, use underscores: `solver_dt`.
"""
function Config(config_path::AbstractString; kwargs...)::Config
    toml = from_toml(Toml, config_path; kwargs...)
    dir = dirname(normpath(config_path))
    Config(toml, dir)
end

Base.getproperty(config::Config, sym::Symbol) = getproperty(getfield(config, :toml), sym)

function Base.show(io::IO, c::Config)
    println(io, "Ribasim Config")
    for field in fieldnames(typeof(c))
        f = getfield(c, field)
        f === nothing || println(io, "\t$field\t= $f")
    end
end

function Base.show(io::IO, c::TableOption)
    first = true
    for field in fieldnames(typeof(c))
        f = getfield(c, field)
        if f !== nothing
            first && (first = false; println(io))
            println(io, "\t\t$field\t= $f")
        end
    end
end

end # module config