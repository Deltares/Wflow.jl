abstract type AbstractRiverFlowModel{T} end
abstract type AbstractOverlandFlowModel{T} end
abstract type AbstractSubsurfaceFlowModel end

abstract type AbstractRoutingMethod end
struct AccucapacityFlux <: AbstractRoutingMethod end
struct KinematicWave <: AbstractRoutingMethod end
struct KinematicWaveStaggered <: AbstractRoutingMethod end
struct LocalInertial <: AbstractRoutingMethod end

struct NoSubsurfaceFlow <: AbstractSubsurfaceFlowModel end
struct NoOverlandFlow{T} <: AbstractOverlandFlowModel{T} end
struct NoRiverFlow{T} <: AbstractRiverFlowModel{T} end

"""
Struct for storing routing model components overland flow `overland_flow`, river flow
`river_flow` and subsurface flow `subsurface_flow`.
"""
@kwdef struct Routing{
    O <: AbstractOverlandFlowModel{<:AbstractRoutingMethod},
    R <: AbstractRiverFlowModel{<:AbstractRoutingMethod},
    S <: AbstractSubsurfaceFlowModel,
}
    overland_flow::O = NoOverlandFlow()
    river_flow::R = NoRiverFlow()
    subsurface_flow::S = NoSubsurfaceFlow()
end
