abstract type AbstractRiverErosionModel end

"Struct for storing river bed and bank erosion model variables"
@with_kw struct RiverErosionModelVariables
    # Potential river bed erosion rate [t dt-1]
    bed::Vector{Float}
    # Potential river bank erosion rate [t dt-1]
    bank::Vector{Float}
end

"Initialize river bed and bank erosion model variables"
function RiverErosionModelVariables(
    n::Int;
    bed::Vector{Float} = fill(MISSING_VALUE, n),
    bank::Vector{Float} = fill(MISSING_VALUE, n),
)
    return RiverErosionModelVariables(; bed = bed, bank = bank)
end

"Struct for storing river erosion model boundary conditions"
@with_kw struct RiverErosionBC
    # Waterlevel [m]
    waterlevel::Vector{Float}
end

"Initialize river erosion model boundary conditions"
function RiverErosionBC(n::Int; waterlevel::Vector{Float} = fill(MISSING_VALUE, n))
    return RiverErosionBC(; waterlevel = waterlevel)
end

"Struct for storing river erosion model parameters"
@with_kw struct RiverErosionParameters
    # Mean diameter [mm] in the river bed/bank
    d50::Vector{Float}
end

"Julian and Torres river erosion model"
@with_kw struct RiverErosionJulianTorresModel <: AbstractRiverErosionModel
    boundary_conditions::RiverErosionBC
    parameters::RiverErosionParameters
    variables::RiverErosionModelVariables
end

"Initialize Julian and Torres river erosion parameters"
function RiverErosionParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "river_bottom-and-bank_sediment__median_diameter")
    d50 = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
    river_parameters = RiverErosionParameters(; d50 = d50)

    return river_parameters
end

"Initialize Julian and Torres river erosion model"
function RiverErosionJulianTorresModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    vars = RiverErosionModelVariables(n)
    params = RiverErosionParameters(dataset, config, indices)
    bc = RiverErosionBC(n)
    model = RiverErosionJulianTorresModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
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
    dt::Float,
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