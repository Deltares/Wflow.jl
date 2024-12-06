abstract type AbstractRiverErosionModel{T} end

"Struct for storing river bed and bank erosion model variables"
@get_units @grid_loc @with_kw struct RiverErosionModelVariables{T}
    # Potential river bed erosion
    bed::Vector{T} | "t dt-1"
    # Potential river bank erosion
    bank::Vector{T} | "t dt-1"
end

"Initialize river bed and bank erosion model variables"
function RiverErosionModelVariables(
    n;
    bed::Vector{T} = fill(mv, n),
    bank::Vector{T} = fill(mv, n),
) where {T}
    return RiverErosionModelVariables{T}(; bed = bed, bank = bank)
end

"Struct for storing river erosion model boundary conditions"
@get_units @grid_loc @with_kw struct RiverErosionBC{T}
    # Waterlevel
    waterlevel::Vector{T} | "t dt-1"
end

"Initialize river erosion model boundary conditions"
function RiverErosionBC(n; waterlevel::Vector{T} = fill(mv, n)) where {T}
    return RiverErosionBC{T}(; waterlevel = waterlevel)
end

"Struct for storing river erosion model parameters"
@get_units @grid_loc @with_kw struct RiverErosionParameters{T}
    # Mean diameter in the river bed/bank
    d50::Vector{T} | "mm"
end

"Julian and Torres river erosion model"
@with_kw struct RiverErosionJulianTorresModel{T} <: AbstractRiverErosionModel{T}
    boundary_conditions::RiverErosionBC{T}
    parameters::RiverErosionParameters{T}
    variables::RiverErosionModelVariables{T}
end

"Initialize Julian and Torres river erosion parameters"
function RiverErosionParameters(dataset, config, indices)
    d50 = ncread(
        dataset,
        config,
        "lateral.river.potential_erosion.parameters.d50";
        sel = indices,
        defaults = 0.1,
        type = Float,
    )
    river_parameters = RiverErosionParameters(; d50 = d50)

    return river_parameters
end

"Initialize Julian and Torres river erosion model"
function RiverErosionJulianTorresModel(dataset, config, indices)
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
function update!(model::RiverErosionJulianTorresModel, geometry::RiverGeometry, dt)
    (; waterlevel) = model.boundary_conditions
    (; d50) = model.parameters
    (; bed, bank) = model.variables

    n = length(waterlevel)
    threaded_foreach(1:n; basesize = 1000) do i
        bed[i], bank[i] = river_erosion_julian_torres(
            waterlevel[i],
            d50[i],
            geometry.width[i],
            geometry.length[i],
            geometry.slope[i],
            dt,
        )
    end
end