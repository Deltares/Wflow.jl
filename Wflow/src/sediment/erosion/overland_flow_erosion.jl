abstract type AbstractOverlandFlowErosionModel end

"Struct for storing overland flow erosion model variables"
@with_kw struct OverlandFlowErosionVariables
    n::Int
    # Total soil erosion rate [t dt-1] from overland flow
    soil_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing overland flow erosion model boundary conditions"
@with_kw struct OverlandFlowErosionBC
    n::Int
    # Overland flow [m3 s-1]
    q::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing ANSWERS overland flow erosion model parameters"
@with_kw struct OverlandFlowErosionAnswersParameters
    # Soil erodibility factor [-]
    usle_k::Vector{Float64}
    # Crop management factor [-]
    usle_c::Vector{Float64}
    # ANSWERS overland flow erosion factor [-]
    answers_overland_flow_factor::Vector{Float64}
end

"Initialize ANSWERS overland flow erosion model parameters"
function OverlandFlowErosionAnswersParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    usle_k = ncread(
        dataset,
        config,
        "soil_erosion__usle_k_factor";
        sel = indices,
        defaults = 0.1,
        type = Float64,
    )
    usle_c = ncread(
        dataset,
        config,
        "soil_erosion__usle_c_factor";
        sel = indices,
        defaults = 0.01,
        type = Float64,
    )
    answers_overland_flow_factor = ncread(
        dataset,
        config,
        "soil_erosion__answers_overland_flow_factor";
        sel = indices,
        defaults = 0.9,
        type = Float64,
    )

    answers_parameters =
        OverlandFlowErosionAnswersParameters(; usle_k, usle_c, answers_overland_flow_factor)
    return answers_parameters
end

"ANSWERS overland flow erosion model"
@with_kw struct OverlandFlowErosionAnswersModel <: AbstractOverlandFlowErosionModel
    n::Int
    boundary_conditions::OverlandFlowErosionBC = OverlandFlowErosionBC(; n)
    parameters::OverlandFlowErosionAnswersParameters
    variables::OverlandFlowErosionVariables = OverlandFlowErosionVariables(; n)
end

"Initialize ANSWERS overland flow erosion model"
function OverlandFlowErosionAnswersModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = OverlandFlowErosionAnswersParameters(dataset, config, indices)
    model = OverlandFlowErosionAnswersModel(; n, parameters)
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
    dt::Float64,
)
    (; q) = model.boundary_conditions
    (; usle_k, usle_c, answers_overland_flow_factor) = model.parameters
    (; soil_erosion_rate) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        soil_erosion_rate[i] = overland_flow_erosion_answers(
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
