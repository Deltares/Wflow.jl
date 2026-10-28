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
struct NoRiverFlowModel <: AbstractRiverFlowModel end

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
    river_flow::R = NoRiverFlowModel()
    subsurface_flow::S = NoSubsurfaceFlowModel()
end
