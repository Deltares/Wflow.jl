"Struct to store atmospheric forcing variables"
@with_kw struct AtmosphericForcing
    n_cells::Int
    # Precipitation [m s⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Potential reference evapotranspiration [m s⁻¹]
    potential_evaporation::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Temperature [K]
    temperature::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct to store hydrological forcing variables"
@with_kw struct HydrologicalForcing
    n_cells::Int
    # Rainfall interception by the vegetation [m s⁻¹]
    interception::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Overland flow depth [m]
    waterlevel_land::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Overland flow discharge [m³ s⁻¹]
    q_land::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # River depth [m]
    waterlevel_river::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # River discharge [m³ s⁻¹]
    q_river::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end
