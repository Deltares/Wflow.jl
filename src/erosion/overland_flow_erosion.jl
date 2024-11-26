abstract type AbstractOverlandFlowErosionModel{T} end

struct NoOverlandFlowErosionModel{T} <: AbstractOverlandFlowErosionModel{T} end

## Overland flow structs and functions
@get_units @grid_loc @with_kw struct OverlandFlowErosionVariables{T}
    # Total soil erosion from overland flow
    amount::Vector{T} | "t dt-1"
end

function OverlandFlowErosionVariables(n; amount::Vector{T} = fill(mv, n)) where {T}
    return OverlandFlowErosionVariables{T}(; amount = amount)
end

@get_units @grid_loc @with_kw struct OverlandFlowErosionBC{T}
    # Overland flow [m3 s-1]
    q::Vector{T}
end

function OverlandFlowErosionBC(n; q::Vector{T} = fill(mv, n)) where {T}
    return OverlandFlowErosionBC{T}(; q = q)
end

# ANSWERS specific structs and functions for rainfall erosion
@get_units @grid_loc @with_kw struct OverlandFlowErosionAnswersParameters{T}
    # Soil erodibility factor
    usle_k::Vector{T} | "-"
    # Crop management factor
    usle_c::Vector{T} | "-"
    # Answers overland flow factor
    answers_k::Vector{T} | "-"
end

function OverlandFlowErosionAnswersParameters(nc, config, inds)
    usle_k = ncread(
        nc,
        config,
        "vertical.overland_flow_erosion.parameters.usle_k";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    usle_c = ncread(
        nc,
        config,
        "vertical.overland_flow_erosion.parameters.usle_c";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    answers_k = ncread(
        nc,
        config,
        "vertical.overland_flow_erosion.parameters.answers_k";
        sel = inds,
        defaults = 0.9,
        type = Float,
    )
    answers_parameters = OverlandFlowErosionAnswersParameters(;
        usle_k = usle_k,
        usle_c = usle_c,
        answers_k = answers_k,
    )
    return answers_parameters
end

@with_kw struct OverlandFlowErosionAnswersModel{T} <: AbstractOverlandFlowErosionModel{T}
    boundary_conditions::OverlandFlowErosionBC{T}
    parameters::OverlandFlowErosionAnswersParameters{T}
    variables::OverlandFlowErosionVariables{T}
end

function OverlandFlowErosionAnswersModel(nc, config, inds)
    n = length(inds)
    vars = OverlandFlowErosionVariables(n)
    params = OverlandFlowErosionAnswersParameters(nc, config, inds)
    bc = OverlandFlowErosionBC(n)
    model = OverlandFlowErosionAnswersModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update_boundary_conditions!(
    model::OverlandFlowErosionAnswersModel,
    hydrometeo_forcing::HydrometeoForcing,
)
    (; q) = model.boundary_conditions
    (; q_land) = hydrometeo_forcing
    @. q = q_land
end

function update_boundary_conditions!(
    model::NoOverlandFlowErosionModel,
    hydrometeo_forcing::HydrometeoForcing,
)
    return nothing
end

function update!(model::OverlandFlowErosionAnswersModel, geometry::LandGeometry, ts)
    (; q) = model.boundary_conditions
    (; usle_k, usle_c, answers_k) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = overland_flow_erosion_answers(
            q[i],
            usle_k[i],
            usle_c[i],
            answers_k[i],
            geometry.slope[i],
            geometry.area[i],
            ts,
        )
    end
end

function update!(model::NoOverlandFlowErosionModel)
    return nothing
end

get_overland_flow_erosion(model::NoOverlandFlowErosionModel) = 0.0
get_overland_flow_erosion(model::OverlandFlowErosionAnswersModel) = model.variables.amount