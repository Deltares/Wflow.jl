abstract type AbstractRiverFlowModel end
abstract type AbstractOverlandFlowModel end
abstract type AbstractSubsurfaceFlowModel end

abstract type AbstractRoutingMethod end
abstract type AbstractStaggeredRoutingMethod <: AbstractRoutingMethod end

struct Manning <: AbstractRoutingMethod end
struct KinematicWave <: AbstractRoutingMethod end
struct ManningStaggered <: AbstractStaggeredRoutingMethod end
struct LocalInertial <: AbstractStaggeredRoutingMethod end

struct NoSubsurfaceFlowModel <: AbstractSubsurfaceFlowModel end
struct NoOverlandFlowModel <: AbstractOverlandFlowModel end

struct NoAllocationRiverModel <: AbstractAllocationModel
    n::Int
end

"""
Struct for storing river boundary conditions when river routing is switched off. Lateral
inflow is stored as it can be useful for an external river routing model. The `reservoir`
boundary condition is set to `nothing`.
"""
@with_kw struct NoRiverFlowBC
    n::Int
    # Lateral inflow [m³ s⁻¹]
    inwater::Vector{Float64} = zeros(n)
    # No reservoir model
    reservoir::Nothing = nothing
end

"""
Struct for storing boundary conditions and allocation when river routing is switched off.
River `allocation` is set to `NoAllocationRiverModel`.
"""
@kwdef struct NoRiverFlowModel <: AbstractRiverFlowModel
    boundary_conditions::NoRiverFlowBC
    allocation::NoAllocationRiverModel
end

"""
Struct for storing routing model components overland flow `overland_flow`, river flow
`river_flow` and subsurface flow `subsurface_flow`.
"""
@kwdef struct Routing{
        O <: AbstractOverlandFlowModel,
        R <: AbstractRiverFlowModel,
        S <: AbstractSubsurfaceFlowModel,
    }
    overland_flow::O = NoOverlandFlowModel()
    river_flow::R
    subsurface_flow::S = NoSubsurfaceFlowModel()
end
