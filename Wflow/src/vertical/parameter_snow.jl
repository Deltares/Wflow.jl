const SnowStateType =
    StateType{(snow_storage = 1, snow_water = 2, snow_melt = 3, runoff = 4)}

# SnowVariables + SnowBC + freeze_rate + melt_rate
@kwdef struct NewSnowCache
    # The number of cells [-]
    n::Int
    # Snow storage [mm]
    snow_storage::Vector{Float64} = zeros(n)
    # Liquid water content in the snow pack [mm]
    snow_water::Vector{Float64} = zeros(n)
    # Snow water equivalent (SWE) [mm]
    swe::Vector{Float64} = fill(MISSING_VALUE, n)
    # Runoff from snowpack [mm Δt⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Snow melt [mm Δt⁻¹]
    snow_melt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Lateral snow (SWE) transport from upstreams cells [mm Δt⁻¹]
    snow_in::Vector{Float64} = zeros(n)
    # Lateral snow (SWE) transport out of a cell [mm Δt⁻¹]
    snow_out::Vector{Float64} = zeros(n)
    # Effective_precip::Vector{Float64}
    effective_precip::Vector{Float64} = fill(MISSING_VALUE, n)
    # Precipitation in the form of snow [mm Δt⁻¹]
    snow_precip::Vector{Float64} = zeros(n)
    # Precipitation in the form of rain [mm Δt⁻¹]
    liquid_precip::Vector{Float64} = zeros(n)
    # The rate at which liquid water freezes [mm Δt⁻¹]
    freeze_rate::Vector{Float64} = zeros(n)
    # The rate at which snow melts [mm Δt⁻¹]
    melt_rate::Vector{Float64} = zeros(n)
end

@kwdef struct NewSnowParameters
    cache::NewSnowCache
    properties::SnowHbvParameters
end

@kwdef struct NewSnowModel <: AbstractSnowModel
    # The snow parameters for all cells
    p::NewSnowParameters
    # An integrator object per cell
    integrators::Vector{WflowIntegrator{SnowStateType}}
end