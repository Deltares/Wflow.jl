abstract type AbstractRiverErosionModel{T} end

"Struct for storing river bed and bank erosion model variables"
@with_kw struct RiverErosionModelVariables{T}
    # Potential river bed erosion rate [t dt-1]
    bed::Vector{T}
    # Potential river bank erosion rate [t dt-1]
    bank::Vector{T}
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
@with_kw struct RiverErosionBC{T}
    # Waterlevel [m]
    waterlevel::Vector{T}
end

"Initialize river erosion model boundary conditions"
function RiverErosionBC(n; waterlevel::Vector{T} = fill(mv, n)) where {T}
    return RiverErosionBC{T}(; waterlevel = waterlevel)
end

"Struct for storing river erosion model parameters"
@with_kw struct RiverErosionParameters{T}
    # Mean diameter [μm] in the river bed/bank
    d50::Vector{T}
end

"Julian and Torres river erosion model"
@with_kw struct RiverErosionJulianTorresModel{T} <: AbstractRiverErosionModel{T}
    boundary_conditions::RiverErosionBC{T}
    parameters::RiverErosionParameters{T}
    variables::RiverErosionModelVariables{T}
end

"Initialize Julian and Torres river erosion parameters"
function RiverErosionParameters(dataset, config, indices)
    lens = lens_input_parameter(config, "river_bottom-and-bank_sediment__d50_diameter")
    d50 = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
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