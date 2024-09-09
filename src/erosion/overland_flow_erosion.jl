abstract type AbstractOverlandFlowErosionModel end

struct NoOverlandFlowErosionModel <: AbstractOverlandFlowErosionModel end

## Overland flow structs and functions
@get_units @with_kw struct OverlandFlowErosionModelVars{T}
    # Total soil erosion from overland flow
    overland_flow_erosion::Vector{T} | "t dt-1"
end

function overland_flow_erosion_model_vars(n)
    vars = OverlandFlowErosionModelVars(; overland_flow_erosion = fill(mv, n))
    return vars
end

@get_units @with_kw struct OverlandFlowErosionBC{T}
    # Overland flow
    overland_flow::Vector{T} | "m dt-1"
end

function overland_flow_erosion_bc(n)
    bc = OverlandFlowErosionBC(; overland_flow = fill(mv, n))
    return bc
end

# ANSWERS specific structs and functions for rainfall erosion
@get_units @with_kw struct OverlandFlowErosionAnswersParameters{T}
    # Soil erodibility factor
    usle_k::Vector{T} | "-"
    # Crop management factor
    usle_c::Vector{T} | "-"
    # Answers overland flow factor
    answers_k::Vector{T} | "-"
    # slope
    slope::Vector{T} | "-"
end

@get_units @with_kw struct OverlandFlowErosionAnswersModel{T} <:
                           AbstractOverlandFlowErosionModel
    boundary_conditions::OverlandFlowErosionBC{T} | "-"
    parameters::OverlandFlowErosionAnswersParameters{T} | "-"
    variables::OverlandFlowErosionModelVars{T} | "-"
end

function initialize_answers_params_overland_flow(nc, config, inds)
    usle_k = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.usle_k";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    usle_c = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.usle_c";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    answers_k = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.answers_k";
        sel = inds,
        defaults = 0.9,
        type = Float,
    )
    slope = ncread(nc, config, "vertical.slope"; sel = inds, defaults = 0.01, type = Float)
    answers_parameters = OverlandFlowErosionAnswersParameters(;
        usle_k = usle_k,
        usle_c = usle_c,
        answers_k = answers_k,
        slope = slope,
    )
    return answers_parameters
end

function initialize_answers_overland_flow_erosion_model(nc, config, inds)
    n = length(inds)
    vars = overland_flow_erosion_model_vars(n)
    params = initialize_answers_params_overland_flow(nc, config, inds)
    bc = overland_flow_erosion_bc(n)
    model = OverlandFlowErosionAnswersModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::OverlandFlowErosionAnswersModel, area, dt)
    (; overland_flow) = model.boundary_conditions
    (; usle_k, usle_c, answers_k, slope) = model.parameters
    (; overland_flow_erosion) = model.variables

    n = length(overland_flow)
    threaded_foreach(1:n; basesize = 1000) do i
        overland_flow_erosion[i] = overland_flow_erosion_answers(
            overland_flow[i],
            usle_k[i],
            usle_c[i],
            answers_k[i],
            slope[i],
            area[i],
            dt,
        )
    end
end

function update!(model::NoOverlandFlowErosionModel)
    return nothing
end

get_overland_flow_erosion(model::NoOverlandFlowErosionModel) = 0.0
get_overland_flow_erosion(model::OverlandFlowErosionAnswersModel) =
    model.variables.overland_flow_erosion