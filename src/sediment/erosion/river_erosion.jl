abstract type AbstractRiverErosionModel{T} end

## Potential direct river erosion structs and functions
@get_units @grid_loc @with_kw struct RiverErosionModelVariables{T}
    # Potential river bed erosion
    bed::Vector{T} | "t dt-1"
    # Potential river bank erosion
    bank::Vector{T} | "t dt-1"
end

function RiverErosionModelVariables(
    n;
    bed::Vector{T} = fill(mv, n),
    bank::Vector{T} = fill(mv, n),
) where {T}
    return RiverErosionModelVariables{T}(; bed = bed, bank = bank)
end

@get_units @grid_loc @with_kw struct RiverErosionBC{T}
    # Waterlevel
    waterlevel::Vector{T} | "t dt-1"
end

function RiverErosionBC(n; waterlevel::Vector{T} = fill(mv, n)) where {T}
    return RiverErosionBC{T}(; waterlevel = waterlevel)
end

# Parameters for the Julian Torres river erosion model
@get_units @grid_loc @with_kw struct RiverErosionParameters{T}
    # Mean diameter in the river bed/bank
    d50::Vector{T} | "mm"
end

@with_kw struct RiverErosionJulianTorresModel{T} <: AbstractRiverErosionModel{T}
    boundary_conditions::RiverErosionBC{T}
    parameters::RiverErosionParameters{T}
    variables::RiverErosionModelVariables{T}
end

function RiverErosionParameters(nc, config, inds)
    d50 = ncread(
        nc,
        config,
        "lateral.river.potential_erosion.parameters.d50";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    river_parameters = RiverErosionParameters(; d50 = d50)

    return river_parameters
end

function RiverErosionJulianTorresModel(nc, config, inds)
    n = length(inds)
    vars = RiverErosionModelVariables(n)
    params = RiverErosionParameters(nc, config, inds)
    bc = RiverErosionBC(n)
    model = RiverErosionJulianTorresModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update_boundary_conditions!(
    model::RiverErosionJulianTorresModel,
    hydrological_forcing::HydrologicalForcing,
)
    (; waterlevel) = model.boundary_conditions
    (; waterlevel_river) = hydrological_forcing
    @. waterlevel = waterlevel_river
end

function update!(model::RiverErosionJulianTorresModel, geometry::RiverParameters, ts)
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
            ts,
        )
    end
end