abstract type AbstractRiverErosionModel end

"Struct for storing river bed and bank erosion model variables"
@with_kw struct RiverErosionModelVariables
    n::Int
    # Potential river bed erosion rate [t dt-1]
    bed::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential river bank erosion rate [t dt-1]
    bank::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing river erosion model boundary conditions"
@with_kw struct RiverErosionBC
    n::Int
    # Waterlevel [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing river erosion model parameters"
@with_kw struct RiverErosionParameters
    # Mean diameter [mm] in the river bed/bank
    d50::Vector{Float64}
end

"Julian and Torres river erosion model"
@with_kw struct RiverErosionJulianTorresModel <: AbstractRiverErosionModel
    n::Int
    boundary_conditions::RiverErosionBC = RiverErosionBC(; n)
    parameters::RiverErosionParameters
    variables::RiverErosionModelVariables = RiverErosionModelVariables(; n)
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
        SoilLoss;
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
    n = length(indices)
    parameters = RiverErosionParameters(dataset, config, indices)
    model = RiverErosionJulianTorresModel(; n, parameters)
    return model
end

"Update river erosion model boundary conditions"
function update_boundary_conditions!(
    model::RiverErosionJulianTorresModel,
    hydrological_forcing::HydrologicalForcing,
)
    (; waterlevel) = model.boundary_conditions
    (; waterlevel_river) = hydrological_forcing
    @. waterlevel = waterlevel_river
end

"Update Julian and Torres river erosion model for a single timestep"
function update!(
    model::RiverErosionJulianTorresModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; waterlevel) = model.boundary_conditions
    (; d50) = model.parameters
    (; bed, bank) = model.variables

    n = length(waterlevel)
    threaded_foreach(1:n; basesize = 1000) do i
        bed[i], bank[i] = river_erosion_julian_torres(
            waterlevel[i],
            d50[i],
            parameters.flow_width[i],
            parameters.flow_length[i],
            parameters.slope[i],
            dt,
        )
    end
end
