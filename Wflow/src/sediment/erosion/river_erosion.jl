abstract type AbstractRiverErosionModel end

"Struct for storing river bed and bank erosion model variables"
@with_kw struct RiverErosionModelVariables
    n_cells::Int
    # Potential river bed erosion rate [t dt-1]
    bed::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Potential river bank erosion rate [t dt-1]
    bank::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing river erosion model boundary conditions"
@with_kw struct RiverErosionBC
    n_cells::Int
    # Waterlevel [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing river erosion model parameters"
@with_kw struct RiverErosionParameters
    # Mean diameter [mm] in the river bed/bank
    d50::Vector{Float64}
end

"Julian and Torres river erosion model"
@with_kw struct RiverErosionJulianTorresModel <: AbstractRiverErosionModel
    n_cells::Int
    boundary_conditions::RiverErosionBC = RiverErosionBC(; n_cells)
    parameters::RiverErosionParameters
    variables::RiverErosionModelVariables = RiverErosionModelVariables(; n_cells)
end

"Initialize Julian and Torres river erosion parameters"
function RiverErosionParameters(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    d50 = ncread(
        dataset,
        config,
        "river_bottom_and_bank_sediment__median_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    river_parameters = RiverErosionParameters(; d50)

    return river_parameters
end

"Initialize Julian and Torres river erosion model"
function RiverErosionJulianTorresModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_cells = length(river_indices_2d)
    parameters = RiverErosionParameters(dataset, config, river_indices_2d)
    river_erosion_model = RiverErosionJulianTorresModel(; n_cells, parameters)
    return river_erosion_model
end

"Update river erosion model boundary conditions"
function update_bc_river_erosion_model!(
    river_erosion_model::RiverErosionJulianTorresModel,
    hydrological_forcing::HydrologicalForcing,
)
    (; waterlevel) = river_erosion_model.boundary_conditions
    (; waterlevel_river) = hydrological_forcing
    @. waterlevel = waterlevel_river
end

"Update Julian and Torres river erosion model for a single timestep"
function update_river_erosion_model!(
    river_erosion_model::RiverErosionJulianTorresModel,
    parameters_river::RiverParameters,
    dt::Float64,
)
    (; boundary_conditions, parameters, variables, n_cells) = river_erosion_model
    (; waterlevel) = boundary_conditions
    (; d50) = parameters
    (; bed, bank) = variables

    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        bed[cell_idx], bank[cell_idx] = river_erosion_julian_torres(
            waterlevel[cell_idx],
            d50[cell_idx],
            parameters_river.flow_width[cell_idx],
            parameters_river.flow_length[cell_idx],
            parameters_river.slope[cell_idx],
            dt,
        )
    end
end
