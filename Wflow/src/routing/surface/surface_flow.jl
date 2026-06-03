abstract type AbstractRiverFlowParameters end
abstract type AbstractRiverFlowVariables end

abstract type AbstractOverlandFlowBC end
abstract type AbstractOverlandFlowParameters end
abstract type AbstractOverlandFlowVariables end

"Struct for storing river flow model boundary conditions"
@with_kw struct RiverFlowBC{R <: Union{ReservoirModel, Nothing}}
    n::Int
    # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    external_inflow::Vector{Float64} = zeros(n)
    # Lateral inflow [m³ s⁻¹]
    inwater::Vector{Float64} = zeros(n)
    # Actual cumulative abstraction from external negative inflow [m³]
    actual_external_abstraction_cumulative::Vector{Float64} = zeros(n)
    # Actual average abstraction from external negative inflow [m³ s⁻¹]
    actual_external_abstraction_average::Vector{Float64} = zeros(n)
    # Abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    abstraction::Vector{Float64} = zeros(n)
    # Exchange with floodplain [m³ s⁻¹]
    floodplain_water_exchange::Vector{Float64} = zeros(n)
    # Reservoir model struct of arrays
    reservoir::R
end

"Struct for storing Manning flow parameters"
@with_kw struct ManningFlowParameters
    slope::Vector{Float64}      # Slope [m m⁻¹]
    mannings_n::Vector{Float64} # Manning's roughness [s m⁻⅓]
    alpha_pow::Float64          # Used in the power part of alpha [-]
    alpha_term::Vector{Float64} # Term used in computation of alpha [-]
    alpha::Vector{Float64}      # Constant in momentum equation A = alpha*Q^beta, based on Manning's equation [s3/5 m1/5]
end

@with_kw struct RiverFlowModel{
    T <: AbstractRoutingMethod,
    F <: Union{AbstractFloodPlainModel, Nothing},
    P <: AbstractRiverFlowParameters,
    V <: AbstractRiverFlowVariables,
    A <: AbstractAllocationModel,
} <: AbstractRiverFlowModel
    routing_method::T
    timestepping::TimeStepping
    boundary_conditions::RiverFlowBC
    parameters::P
    variables::V
    floodplain::F
    allocation::A
end

@with_kw struct OverlandFlowModel{
    T <: AbstractRoutingMethod,
    B <: AbstractOverlandFlowBC,
    P <: Union{ManningFlowParameters, AbstractOverlandFlowParameters},
    V <: AbstractOverlandFlowVariables,
} <: AbstractOverlandFlowModel
    routing_method::T
    timestepping::TimeStepping
    boundary_conditions::B
    parameters::P
    variables::V
end

function wetted_perimeter_channel(h::Float64, flow_width::Float64)
    channel_perimeter = 2.0 * h + flow_width
    return channel_perimeter
end
