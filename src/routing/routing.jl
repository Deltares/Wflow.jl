abstract type AbstractSubsurfaceFlowModel end
abstract type AbstractOverlandFlowModel end
abstract type AbstractRiverFlowModel end

struct NoSubsurfaceFlow <: AbstractSubsurfaceFlowModel end
struct NoOverlandFlow <: AbstractOverlandFlowModel end
struct NoRiverFlow <: AbstractRiverFlowModel end

""" 
Struct for storing routing model components overland flow `overland_flow`, river flow
`river_flow` and subsurface flow `subsurface_flow`.
"""
@kwdef struct Routing{
    O <: AbstractOverlandFlowModel,
    R <: AbstractRiverFlowModel,
    S <: AbstractSubsurfaceFlowModel,
}
    overland_flow::O = NoSubsurfaceFlow()
    river_flow::R = NoOverlandFlow()
    subsurface_flow::S = NoSubsurfaceFlow()
end