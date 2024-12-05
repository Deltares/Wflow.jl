abstract type AbstractOverlandFlowErosionModel{T} end

"Struct for storing overland flow erosion model variables"
@get_units @grid_loc @with_kw struct OverlandFlowErosionVariables{T}
    # Total soil erosion from overland flow
    amount::Vector{T} | "t dt-1"
end

"Initialize overland flow erosion model variables"
function OverlandFlowErosionVariables(n; amount::Vector{T} = fill(mv, n)) where {T}
    return OverlandFlowErosionVariables{T}(; amount = amount)
end

"Struct for storing overland flow erosion model boundary conditions"
@get_units @grid_loc @with_kw struct OverlandFlowErosionBC{T}
    # Overland flow [m3 s-1]
    q::Vector{T}
end

"Initialize overland flow erosion model boundary conditions"
function OverlandFlowErosionBC(n; q::Vector{T} = fill(mv, n)) where {T}
    return OverlandFlowErosionBC{T}(; q = q)
end

"Struct for storing ANSWERS overland flow erosion model parameters"
@get_units @grid_loc @with_kw struct OverlandFlowErosionAnswersParameters{T}
    # Soil erodibility factor
    usle_k::Vector{T} | "-"
    # Crop management factor
    usle_c::Vector{T} | "-"
    # Answers overland flow factor
    answers_k::Vector{T} | "-"
end

"Initialize ANSWERS overland flow erosion model parameters"
function OverlandFlowErosionAnswersParameters(dataset, config, indices)
    usle_k = ncread(
        dataset,
        config,
        "vertical.overland_flow_erosion.parameters.usle_k";
        sel = indices,
        defaults = 0.1,
        type = Float,
    )
    usle_c = ncread(
        dataset,
        config,
        "vertical.overland_flow_erosion.parameters.usle_c";
        sel = indices,
        defaults = 0.01,
        type = Float,
    )
    answers_k = ncread(
        dataset,
        config,
        "vertical.overland_flow_erosion.parameters.answers_k";
        sel = indices,
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

"ANSWERS overland flow erosion model"
@with_kw struct OverlandFlowErosionAnswersModel{T} <: AbstractOverlandFlowErosionModel{T}
    boundary_conditions::OverlandFlowErosionBC{T}
    parameters::OverlandFlowErosionAnswersParameters{T}
    variables::OverlandFlowErosionVariables{T}
end

"Initialize ANSWERS overland flow erosion model"
function OverlandFlowErosionAnswersModel(dataset, config, indices)
    n = length(indices)
    vars = OverlandFlowErosionVariables(n)
    params = OverlandFlowErosionAnswersParameters(dataset, config, indices)
    bc = OverlandFlowErosionBC(n)
    model = OverlandFlowErosionAnswersModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update boundary conditions for ANSWERS overland flow erosion model"
function update_boundary_conditions!(
    model::OverlandFlowErosionAnswersModel,
    hydrological_forcing::HydrologicalForcing,
)
    (; q) = model.boundary_conditions
    (; q_land) = hydrological_forcing
    @. q = q_land
end

"Update ANSWERS overland flow erosion model for a single timestep"
function update!(model::OverlandFlowErosionAnswersModel, geometry::LandParameters, dt)
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
            dt,
        )
    end
end