"""
A WflowIntegrator object holds all parameters and variables
needed for numerical integration of an ODE.
"""
@kwdef mutable struct WflowIntegrator{V <: ComponentVector}
    # The fraction of the global time step that has been simulated
    progress::Float64 = 0.0
    # The sub time step
    dt_sub::Float64 = 0.0
    # The state vector
    const u::V
    # The previous value of the state vector
    const uprev::V
    # The time derivative of the state vector,
    # generally instantaneous fluxes
    const du::V
end

function WflowIntegrator(u0::V) where {V <: ComponentVector}
    return WflowIntegrator(; u = u0, uprev = copy(u0), du = zero(u0))
end

const StateType{A} = ComponentVector{Float64, Vector{Float64}, Tuple{Axis{A}}}
const SnowStateType = StateType{(
    snow_storage = 1,
    snow_water = 2,
    cumulative_snow_melt = 3,
    cumulative_runoff = 4,
)}

@kwdef struct NewSnowCache
    # The number of cells [-]
    n::Int
    # Snow water equivalent (SWE) [mm]
    swe::Vector{Float64} = fill(MISSING_VALUE, n)
    # Snow melt [mm Δt⁻¹]
    snow_melt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Runoff from snowpack [mm Δt⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Lateral snow (SWE) transport from upstreams cells [mm Δt⁻¹]
    snow_in::Vector{Float64} = zeros(n)
    # Lateral snow (SWE) transport out of a cell [mm Δt⁻¹]
    snow_out::Vector{Float64} = zeros(n)
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