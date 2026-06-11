abstract type AbstractRiverErosionModel end

"Struct for storing river bed and bank erosion model variables"
@with_kw struct RiverErosionModelVariables
    n_river::Int
    # Potential river bed erosion rate [kg s⁻¹]
    bed::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Potential river bank erosion rate [kg s⁻¹]
    bank::Vector{Float64} = fill(MISSING_VALUE, n_river)
end

"Struct for storing river erosion model boundary conditions"
@with_kw struct RiverErosionBC
    n_river::Int
    # Waterlevel [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n_river)
end

"Struct for storing river erosion model parameters"
@with_kw struct RiverErosionParameters
    # Mean diameter [m] in the river bed/bank
    d50::Vector{Float64}
end

"Julian and Torres river erosion model"
@with_kw struct RiverErosionJulianTorresModel <: AbstractRiverErosionModel
    n_river::Int
    boundary_conditions::RiverErosionBC = RiverErosionBC(; n_river)
    parameters::RiverErosionParameters
    variables::RiverErosionModelVariables = RiverErosionModelVariables(; n_river)
end

"Initialize Julian and Torres river erosion parameters"
function RiverErosionParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    d50 = ncread(
        dataset,
        config,
        "river_bottom_and_bank_sediment__median_diameter",
        SoilLossModel;
        sel = indices,
    )
    river_parameters = RiverErosionParameters(; d50)

    return river_parameters
end

"Initialize Julian and Torres river erosion model"
function RiverErosionJulianTorresModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n_river = length(indices)
    parameters = RiverErosionParameters(dataset, config, indices)
    river_erosion_model = RiverErosionJulianTorresModel(; n_river, parameters)
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
    (; boundary_conditions, parameters, variables) = river_erosion_model
    (; waterlevel) = boundary_conditions
    (; d50) = parameters
    (; bed, bank) = variables

    threaded_foreach(eachindex(waterlevel); basesize = 1000) do river_idx
        bed[river_idx], bank[river_idx] = river_erosion_julian_torres(
            waterlevel[river_idx],
            d50[river_idx],
            parameters_river.flow_width[river_idx],
            parameters_river.flow_length[river_idx],
            parameters_river.slope[river_idx],
            dt,
        )
    end
end
