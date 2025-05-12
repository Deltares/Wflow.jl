abstract type AbstractOverlandFlowErosionModel end

"Struct for storing overland flow erosion model variables"
@with_kw struct OverlandFlowErosionVariables
    # Total soil erosion rate [t dt-1] from overland flow
    amount::Vector{Float}
end

"Initialize overland flow erosion model variables"
function OverlandFlowErosionVariables(
    n::Int;
    amount::Vector{Float} = fill(MISSING_VALUE, n),
)
    return OverlandFlowErosionVariables(; amount = amount)
end

"Struct for storing overland flow erosion model boundary conditions"
@with_kw struct OverlandFlowErosionBC
    # Overland flow [m3 s-1]
    q::Vector{Float}
end

"Initialize overland flow erosion model boundary conditions"
function OverlandFlowErosionBC(n::Int; q::Vector{Float} = fill(MISSING_VALUE, n))
    return OverlandFlowErosionBC(; q = q)
end

"Struct for storing ANSWERS overland flow erosion model parameters"
@with_kw struct OverlandFlowErosionAnswersParameters
    # Soil erodibility factor [-]
    usle_k::Vector{Float}
    # Crop management factor [-]
    usle_c::Vector{Float}
    # ANSWERS overland flow erosion factor [-]
    answers_overland_flow_factor::Vector{Float}
end

"Initialize ANSWERS overland flow erosion model parameters"
function OverlandFlowErosionAnswersParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "soil_erosion__usle_k_factor")
    usle_k = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
    lens = lens_input_parameter(config, "soil_erosion__usle_c_factor")
    usle_c = ncread(dataset, config, lens; sel = indices, defaults = 0.01, type = Float)
    lens = lens_input_parameter(config, "soil_erosion__answers_overland_flow_factor")
    answers_overland_flow_factor =
        ncread(dataset, config, lens; sel = indices, defaults = 0.9, type = Float)

    answers_parameters = OverlandFlowErosionAnswersParameters(;
        usle_k = usle_k,
        usle_c = usle_c,
        answers_overland_flow_factor = answers_overland_flow_factor,
    )
    return answers_parameters
end

"ANSWERS overland flow erosion model"
@with_kw struct OverlandFlowErosionAnswersModel <: AbstractOverlandFlowErosionModel
    boundary_conditions::OverlandFlowErosionBC
    parameters::OverlandFlowErosionAnswersParameters
    variables::OverlandFlowErosionVariables
end

"Initialize ANSWERS overland flow erosion model"
function OverlandFlowErosionAnswersModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
function update!(
    model::OverlandFlowErosionAnswersModel,
    geometry::LandParameters,
    dt::Float,
)
    (; q) = model.boundary_conditions
    (; usle_k, usle_c, answers_overland_flow_factor) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = overland_flow_erosion_answers(
            q[i],
            usle_k[i],
            usle_c[i],
            answers_overland_flow_factor[i],
            geometry.slope[i],
            geometry.area[i],
            dt,
        )
    end
end