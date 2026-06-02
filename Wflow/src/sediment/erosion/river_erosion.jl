abstract type AbstractRiverErosionModel end

"Struct for storing river bed and bank erosion model variables"
@kwdef struct RiverErosionModelVariables
    n::Int
    # Potential river bed erosion rate [kg s⁻¹]
    bed::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential river bank erosion rate [kg s⁻¹]
    bank::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing river erosion model boundary conditions"
@kwdef struct RiverErosionBC
    n::Int
    # Waterlevel [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing river erosion model parameters"
@with_data_lookup struct RiverErosionParameters
    # Mean diameter [m] in the river bed/bank
    "river_bottom_and_bank_sediment__median_diameter"
    d50::Vector{Float64}
end

"Julian and Torres river erosion model"
@kwdef struct RiverErosionJulianTorresModel <: AbstractRiverErosionModel
    n::Int
    boundary_conditions::RiverErosionBC = RiverErosionBC(; n)
    parameters::RiverErosionParameters
    variables::RiverErosionModelVariables = RiverErosionModelVariables(; n)
end

"Initialize Julian and Torres river erosion parameters"
function RiverErosionParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    d50 = ncread(
        dataset,
        config,
        "river_bottom_and_bank_sediment__median_diameter",
        SoilLossModel;
        sel = indices,
    )
    river_parameters = RiverErosionParameters(data_lookup; d50)

    return river_parameters
end

"Initialize Julian and Torres river erosion model"
function RiverErosionJulianTorresModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    n = length(indices)
    parameters = RiverErosionParameters(dataset, config, indices; data_lookup)
    river_erosion_model = RiverErosionJulianTorresModel(; n, parameters)
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

    threaded_foreach(eachindex(waterlevel); basesize = 1000) do i
        bed[i], bank[i] = river_erosion_julian_torres(
            waterlevel[i],
            d50[i],
            parameters_river.flow_width[i],
            parameters_river.flow_length[i],
            parameters_river.slope[i],
            dt,
        )
    end
end
