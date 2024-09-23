abstract type AbstractOverlandFlowErosionModel end

struct NoOverlandFlowErosionModel <: AbstractOverlandFlowErosionModel end

## Overland flow structs and functions
@get_units @with_kw struct OverlandFlowErosionModelVars{T}
    # Total soil erosion from overland flow
    amount::Vector{T} | "t dt-1"
end

function overland_flow_erosion_model_vars(n)
    vars = OverlandFlowErosionModelVars(; amount = fill(mv, n))
    return vars
end

# ANSWERS specific structs and functions for rainfall erosion
@get_units @with_kw struct OverlandFlowErosionAnswersParameters{T}
    # Soil erodibility factor
    usle_k::Vector{T} | "-"
    # Crop management factor
    usle_c::Vector{T} | "-"
    # Answers overland flow factor
    answers_k::Vector{T} | "-"
end

@get_units @with_kw struct OverlandFlowErosionAnswersModel{T} <:
                           AbstractOverlandFlowErosionModel
    parameters::OverlandFlowErosionAnswersParameters{T} | "-"
    variables::OverlandFlowErosionModelVars{T} | "-"
end

function initialize_answers_params_overland_flow(nc, config, inds)
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

function initialize_answers_overland_flow_erosion_model(nc, config, inds)
    n = length(inds)
    vars = overland_flow_erosion_model_vars(n)
    params = initialize_answers_params_overland_flow(nc, config, inds)
    model = OverlandFlowErosionAnswersModel(; parameters = params, variables = vars)
    return model
end

function update!(
    model::OverlandFlowErosionAnswersModel,
    hydrometeo_forcing::HydrometeoForcing,
    geometry::LandGeometry,
    ts,
)
    (; q_land) = hydrometeo_forcing
    (; usle_k, usle_c, answers_k) = model.parameters
    (; amount) = model.variables

    n = length(q_land)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = overland_flow_erosion_answers(
            q_land[i],
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