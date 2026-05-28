abstract type AbstractOverlandFlowErosionModel end

"Struct for storing overland flow erosion model variables"
@with_kw struct OverlandFlowErosionVariables
    n_cells::Int
    # Total soil erosion rate [t dt-1] from overland flow
    soil_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing overland flow erosion model boundary conditions"
@with_kw struct OverlandFlowErosionBC
    n_cells::Int
    # Overland flow [m3 s-1]
    q::Vector{Float64} = fill(MISSING_VALUE, n_cells)
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
    land_indices_2d::Vector{CartesianIndex{2}},
)
    usle_k = ncread(
        dataset,
        config,
        "soil_erosion__usle_k_factor",
        SoilLossModel;
        sel = land_indices_2d,
    )
    usle_c = ncread(
        dataset,
        config,
        "soil_erosion__usle_c_factor",
        SoilLossModel;
        sel = land_indices_2d,
    )
    answers_overland_flow_factor = ncread(
        dataset,
        config,
        "soil_erosion__answers_overland_flow_factor",
        SoilLossModel;
        sel = land_indices_2d,
    )

    answers_parameters =
        OverlandFlowErosionAnswersParameters(; usle_k, usle_c, answers_overland_flow_factor)
    return answers_parameters
end

"ANSWERS overland flow erosion model"
@with_kw struct OverlandFlowErosionAnswersModel <: AbstractOverlandFlowErosionModel
    n_cells::Int
    boundary_conditions::OverlandFlowErosionBC = OverlandFlowErosionBC(; n_cells)
    parameters::OverlandFlowErosionAnswersParameters
    variables::OverlandFlowErosionVariables = OverlandFlowErosionVariables(; n_cells)
end

"Initialize ANSWERS overland flow erosion model"
function OverlandFlowErosionAnswersModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
)
    n_cells = length(land_indices_2d)
    parameters = OverlandFlowErosionAnswersParameters(dataset, config, land_indices_2d)
    overland_flow_erosion_model = OverlandFlowErosionAnswersModel(; n_cells, parameters)
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
    (; n_cells) = overland_flow_erosion_model
    (; q) = overland_flow_erosion_model.boundary_conditions
    (; usle_k, usle_c, answers_overland_flow_factor) =
        overland_flow_erosion_model.parameters
    (; soil_erosion_rate) = overland_flow_erosion_model.variables

    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        soil_erosion_rate[cell_idx] = overland_flow_erosion_answers(
            q[cell_idx],
            usle_k[cell_idx],
            usle_c[cell_idx],
            answers_overland_flow_factor[cell_idx],
            geometry.slope[cell_idx],
            geometry.area[cell_idx],
            dt,
        )
    end
end
