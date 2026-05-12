abstract type AbstractOverlandFlowErosionModel end

"Struct for storing overland flow erosion model variables"
@with_data_lookup struct OverlandFlowErosionVariables
    n::Int
    # Total soil erosion rate [t dt-1] from overland flow
    "overland_flow_soil_erosion__mass_flow_rate"
    soil_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing overland flow erosion model boundary conditions"
@with_kw struct OverlandFlowErosionBC
    n::Int
    # Overland flow [m3 s-1]
    q::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing ANSWERS overland flow erosion model parameters"
@with_data_lookup struct OverlandFlowErosionAnswersParameters
    # Soil erodibility factor [-]
    "soil_erosion__usle_k_factor"
    usle_k::Vector{Float64}
    # Crop management factor [-]
    "soil_erosion__usle_c_factor"
    usle_c::Vector{Float64}
    # ANSWERS overland flow erosion factor [-]
    "soil_erosion__answers_overland_flow_factor"
    answers_overland_flow_factor::Vector{Float64}
end

"Initialize ANSWERS overland flow erosion model parameters"
function OverlandFlowErosionAnswersParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    usle_k =
        ncread(dataset, config, "soil_erosion__usle_k_factor", SoilLossModel; sel = indices)
    usle_c =
        ncread(dataset, config, "soil_erosion__usle_c_factor", SoilLossModel; sel = indices)
    answers_overland_flow_factor = ncread(
        dataset,
        config,
        "soil_erosion__answers_overland_flow_factor",
        SoilLossModel;
        sel = indices,
    )

    answers_parameters = OverlandFlowErosionAnswersParameters(
        data_lookup;
        usle_k,
        usle_c,
        answers_overland_flow_factor,
    )
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
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    n = length(indices)
    parameters = OverlandFlowErosionAnswersParameters(dataset, config, indices; data_lookup)
    overland_flow_erosion_model = OverlandFlowErosionAnswersModel(; n, parameters)
    return overland_flow_erosion_model
end

"Update boundary conditions for ANSWERS overland flow erosion model"
function update_bc_overland_flow_erosion_model!(
    overland_flow_erosion_model::OverlandFlowErosionAnswersModel,
    hydrological_forcing::HydrologicalForcing,
)
    (; q) = overland_flow_erosion_model.boundary_conditions
    (; q_land) = hydrological_forcing
    @. q = q_land
end

"Update ANSWERS overland flow erosion model for a single timestep"
function update_overland_flow_erosion_model!(
    overland_flow_erosion_model::OverlandFlowErosionAnswersModel,
    geometry::LandParameters,
    dt::Float64,
)
    (; q) = overland_flow_erosion_model.boundary_conditions
    (; usle_k, usle_c, answers_overland_flow_factor) =
        overland_flow_erosion_model.parameters
    (; soil_erosion_rate) = overland_flow_erosion_model.variables

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
