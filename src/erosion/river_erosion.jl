abstract type AbstractRiverErosionModel end

## Potential direct river erosion structs and functions
@get_units @with_kw struct RiverErosionModelVars{T}
    # Potential river bed erosion
    bed::Vector{T} | "t dt-1"
    # Potential river bank erosion
    bank::Vector{T} | "t dt-1"
end

function river_erosion_model_vars(n)
    vars = RiverErosionModelVars(; bed = fill(mv, n), bank = fill(mv, n))
    return vars
end

@get_units @with_kw struct RiverErosionBC{T}
    # Waterlevel
    waterlevel::Vector{T} | "t dt-1"
end

function river_erosion_bc(n)
    bc = RiverErosionBC(; waterlevel = fill(mv, n))
    return bc
end

# Parameters for the Julian Torres river erosion model
@get_units @with_kw struct RiverErosionParameters{T}
    # Mean diameter in the river bed/bank
    d50::Vector{T} | "mm"
end

@get_units @with_kw struct RiverErosionJulianTorresModel{T} <: AbstractRiverErosionModel
    boundary_conditions::RiverErosionBC{T} | "-"
    parameters::RiverErosionParameters{T} | "-"
    variables::RiverErosionModelVars{T} | "-"
end

function initialize_river_erosion_params(nc, config, inds)
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

function initialize_river_erosion_julian_torres_model(nc, config, inds)
    n = length(inds)
    vars = river_erosion_model_vars(n)
    params = initialize_river_erosion_params(nc, config, inds)
    bc = river_erosion_bc(n)
    model = RiverErosionJulianTorresModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::RiverErosionJulianTorresModel, geometry::RiverGeometry, ts)
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