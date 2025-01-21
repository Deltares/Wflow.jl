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

function Routing(
    config::Config,
    dataset::NCDataset,
    data::NamedTuple,
    todo::NamedTuple,
    type::AbstractSbmModelType,
)::Tuple{Routing, NamedTuple}
    (; sub_catchment_data, cell_data) = data
    (; ldd) = sub_catchment_data
    (; x_length, y_length) = cell_data
    flow_length = map(get_flow_length, ldd, x_length, y_length)
    flow_width = (x_length .* y_length) ./ flow_length
    flow_data = (; flow_length, flow_width)
    data = (; data..., flow_data)

    subsurface_flow = get_subsurface_flow(config, dataset, data, todo, type)
    overland_flow = get_overland_flow(config, dataset, data, type)
    river_flow = get_river_flow(config, dataset, data, type)

    routing = Routing(; subsurface_flow, overland_flow, river_flow)
    data = (; data..., routing)

    return routing, data
end