abstract type AbstractTransportCapacityModel{T} end

"Struct to store total transport capacity model variables"
@get_units @grid_loc @with_kw struct TransportCapacityModelVariables{T}
    # Total sediment transport capacity
    amount::Vector{T} | "t dt-1"
end

"Initialize total transport capacity model variables"
function TransportCapacityModelVariables(n; amount::Vector{T} = fill(mv, n)) where {T}
    return TransportCapacityModelVariables{T}(; amount = amount)
end

"Struct to store total transport capacity model boundary conditions"
@get_units @grid_loc @with_kw struct TransportCapacityBC{T}
    # Discharge
    q::Vector{T} | "m3 s-1"
    # Flow depth
    waterlevel::Vector{T} | "m"
end

"Initialize total transport capacity model boundary conditions"
function TransportCapacityBC(
    n;
    q::Vector{T} = fill(mv, n),
    waterlevel::Vector{T} = fill(mv, n),
) where {T}
    return TransportCapacityBC{T}(; q = q, waterlevel = waterlevel)
end

"Update total transport capacity model boundary conditions"
function update_boundary_conditions!(
    model::AbstractTransportCapacityModel,
    hydrological_forcing::HydrologicalForcing,
    model_type::Symbol,
)
    (; q, waterlevel) = model.boundary_conditions
    (; q_land, waterlevel_land, q_river, waterlevel_river) = hydrological_forcing

    if model_type == :land
        @. q = q_land
        @. waterlevel = waterlevel_land
    elseif model_type == :river
        @. q = q_river
        @. waterlevel = waterlevel_river
    end
end

##################### Overland Flow #####################

"Struct to store Govers overland flow transport capacity model parameters"
@get_units @grid_loc @with_kw struct TransportCapacityGoversParameters{T}
    # Particle density
    density::Vector{T} | "kg m-3"
    # Govers transport capacity coefficient
    c_govers::Vector{T} | "-"
    # Govers transport capacity exponent
    n_govers::Vector{T} | "-"
end

"Initialize Govers overland flow transport capacity model parameters"
function TransportCapacityGoversParameters(dataset, config, indices)
    lens = lens_input_parameter(config, "land_surface_sediment__particle_density")
    density = ncread(dataset, config, lens; sel = indices, defaults = 2650.0, type = Float)
    lens = lens_input_parameter(
        config,
        "land_surface_water_sediment__govers_transport_capacity_coefficient",
    )
    c_govers =
        ncread(dataset, config, lens; sel = indices, defaults = 0.000505, type = Float)
    lens = lens_input_parameter(
        config,
        "land_surface_water_sediment__govers_transport_capacity_exponent",
    )
    n_govers = ncread(dataset, config, lens; sel = indices, defaults = 4.27, type = Float)
    tc_parameters = TransportCapacityGoversParameters(;
        density = density,
        c_govers = c_govers,
        n_govers = n_govers,
    )

    return tc_parameters
end

"Govers overland flow transport capacity model"
@with_kw struct TransportCapacityGoversModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityGoversParameters{T}
    variables::TransportCapacityModelVariables{T}
end

"Initialize Govers overland flow transport capacity model"
function TransportCapacityGoversModel(dataset, config, indices)
    n = length(indices)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityGoversParameters(dataset, config, indices)
    bc = TransportCapacityBC(n)
    model = TransportCapacityGoversModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update Govers overland flow transport capacity model for a single timestep"
function update!(
    model::TransportCapacityGoversModel,
    geometry::LandGeometry,
    waterbodies,
    rivers,
    dt,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, c_govers, n_govers) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_govers(
            q[i],
            waterlevel[i],
            c_govers[i],
            n_govers[i],
            density[i],
            geometry.slope[i],
            geometry.width[i],
            waterbodies[i],
            rivers[i],
            dt,
        )
    end
end

"Struct to store Yalin overland flow transport capacity model parameters"
@get_units @grid_loc @with_kw struct TransportCapacityYalinParameters{T}
    # Particle density
    density::Vector{T} | "kg m-3"
    # Particle mean diameter
    d50::Vector{T} | "mm"
end

"Initialize Yalin overland flow transport capacity model parameters"
function TransportCapacityYalinParameters(dataset, config, indices)
    lens = lens_input_parameter(config, "land_surface_sediment__particle_density")
    density = ncread(dataset, config, lens; sel = indices, defaults = 2650.0, type = Float)

    lens = lens_input_parameter(config, "land_surface_sediment__d50_diameter")
    d50 = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)

    tc_parameters = TransportCapacityYalinParameters(; density = density, d50 = d50)

    return tc_parameters
end

"Yalin overland flow transport capacity model"
@with_kw struct TransportCapacityYalinModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityYalinParameters{T}
    variables::TransportCapacityModelVariables{T}
end

"Initialize Yalin overland flow transport capacity model"
function TransportCapacityYalinModel(dataset, config, indices)
    n = length(indices)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityYalinParameters(dataset, config, indices)
    bc = TransportCapacityBC(n)
    model = TransportCapacityYalinModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update Yalin overland flow transport capacity model for a single timestep"
function update!(
    model::TransportCapacityYalinModel,
    geometry::LandGeometry,
    waterbodies,
    rivers,
    dt,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_yalin(
            q[i],
            waterlevel[i],
            density[i],
            d50[i],
            geometry.slope[i],
            geometry.width[i],
            waterbodies[i],
            rivers[i],
            dt,
        )
    end
end

"Struct to store Yalin differentiated overland flow transport capacity model variables"
@get_units @grid_loc @with_kw struct TransportCapacityYalinDifferentiationModelVariables{T}
    # Total sediment transport capacity
    amount::Vector{T} | "t dt-1"
    # Transport capacity clay
    clay::Vector{T} | "t dt-1"
    # Transport capacity silt
    silt::Vector{T} | "t dt-1"
    # Transport capacity sand
    sand::Vector{T} | "t dt-1"
    # Transport capacity small aggregates
    sagg::Vector{T} | "t dt-1"
    # Transport capacity large aggregates
    lagg::Vector{T} | "t dt-1"
end

"Initialize Yalin differentiated overland flow transport capacity model variables"
function TransportCapacityYalinDifferentiationModelVariables(
    n;
    amount::Vector{T} = fill(mv, n),
    clay::Vector{T} = fill(mv, n),
    silt::Vector{T} = fill(mv, n),
    sand::Vector{T} = fill(mv, n),
    sagg::Vector{T} = fill(mv, n),
    lagg::Vector{T} = fill(mv, n),
) where {T}
    return TransportCapacityYalinDifferentiationModelVariables{T}(;
        amount = amount,
        clay = clay,
        silt = silt,
        sand = sand,
        sagg = sagg,
        lagg = lagg,
    )
end

"Struct to store Yalin differentiated overland flow transport capacity model parameters"
@get_units @grid_loc @with_kw struct TransportCapacityYalinDifferentiationParameters{T}
    # Particle density
    density::Vector{T} | "kg m-3"
    # Clay mean diameter
    dm_clay::Vector{T} | "µm"
    # Silt mean diameter
    dm_silt::Vector{T} | "µm"
    # Sand mean diameter
    dm_sand::Vector{T} | "µm"
    # Small aggregates mean diameter
    dm_sagg::Vector{T} | "µm"
    # Large aggregates mean diameter
    dm_lagg::Vector{T} | "µm"
end

"Initialize Yalin differentiated overland flow transport capacity model parameters"
function TransportCapacityYalinDifferentiationParameters(dataset, config, indices)
    lens = lens_input_parameter(config, "land_surface_sediment__particle_density")
    density = ncread(dataset, config, lens; sel = indices, defaults = 2650.0, type = Float)
    lens = lens_input_parameter(config, "clay__d50_diameter")
    dm_clay = ncread(dataset, config, lens; sel = indices, defaults = 2.0, type = Float)
    lens = lens_input_parameter(config, "silt__d50_diameter")
    dm_silt = ncread(dataset, config, lens; sel = indices, defaults = 10.0, type = Float)
    lens = lens_input_parameter(config, "sand__d50_diameter")
    dm_sand = ncread(dataset, config, lens; sel = indices, defaults = 200.0, type = Float)
    lens = lens_input_parameter(config, "sediment_aggregates~small__d50_diameter")
    dm_sagg = ncread(dataset, config, lens; sel = indices, defaults = 30.0, type = Float)
    lens = lens_input_parameter(config, "sediment_aggregates~large__d50_diameter")
    dm_lagg = ncread(dataset, config, lens; sel = indices, defaults = 500.0, type = Float)

    tc_parameters = TransportCapacityYalinDifferentiationParameters(;
        density = density,
        dm_clay = dm_clay,
        dm_silt = dm_silt,
        dm_sand = dm_sand,
        dm_sagg = dm_sagg,
        dm_lagg = dm_lagg,
    )

    return tc_parameters
end

"Yalin differentiated overland flow transport capacity model"
@with_kw struct TransportCapacityYalinDifferentiationModel{T} <:
                AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityYalinDifferentiationParameters{T}
    variables::TransportCapacityYalinDifferentiationModelVariables{T}
end

"Initialize Yalin differentiated overland flow transport capacity model"
function TransportCapacityYalinDifferentiationModel(dataset, config, indices)
    n = length(indices)
    vars = TransportCapacityYalinDifferentiationModelVariables(n)
    params = TransportCapacityYalinDifferentiationParameters(dataset, config, indices)
    bc = TransportCapacityBC(n)
    model = TransportCapacityYalinDifferentiationModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update Yalin differentiated overland flow transport capacity model for a single timestep"
function update!(
    model::TransportCapacityYalinDifferentiationModel,
    geometry::LandGeometry,
    waterbodies,
    rivers,
    dt,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg) = model.parameters
    (; amount, clay, silt, sand, sagg, lagg) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        dtot = transportability_yalin_differentiation(
            waterlevel[i],
            density[i],
            dm_clay[i],
            dm_silt[i],
            dm_sand[i],
            dm_sagg[i],
            dm_lagg[i],
            geometry.slope[i],
        )
        clay[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_clay[i],
            geometry.slope[i],
            geometry.width[i],
            waterbodies[i],
            rivers[i],
            dtot,
            dt,
        )
        silt[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_silt[i],
            geometry.slope[i],
            geometry.width[i],
            waterbodies[i],
            rivers[i],
            dtot,
            dt,
        )
        sand[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_sand[i],
            geometry.slope[i],
            geometry.width[i],
            waterbodies[i],
            rivers[i],
            dtot,
            dt,
        )
        sagg[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_sagg[i],
            geometry.slope[i],
            geometry.width[i],
            waterbodies[i],
            rivers[i],
            dtot,
            dt,
        )
        lagg[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_lagg[i],
            geometry.slope[i],
            geometry.width[i],
            waterbodies[i],
            rivers[i],
            dtot,
            dt,
        )
        amount[i] = clay[i] + silt[i] + sand[i] + sagg[i] + lagg[i]
    end
end

"Struct to store common river transport capacity model parameters"
@get_units @grid_loc @with_kw struct TransportCapacityRiverParameters{T}
    # Particle density
    density::Vector{T} | "kg m-3"
    # Particle mean diameter
    d50::Vector{T} | "mm"
end

"Initialize common river transport capacity model parameters"
function TransportCapacityRiverParameters(dataset, config, indices)
    lens = lens_input_parameter(config, "river_sediment__particle_density")
    density = ncread(dataset, config, lens; sel = indices, defaults = 2650.0, type = Float)
    lens = lens_input_parameter(config, "river_sediment__d50_diameter")
    d50 = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)

    tc_parameters = TransportCapacityRiverParameters(; density = density, d50 = d50)

    return tc_parameters
end

"Struct to store Bagnold transport capacity model parameters"
@get_units @grid_loc @with_kw struct TransportCapacityBagnoldParameters{T}
    # Bagnold transport capacity coefficient
    c_bagnold::Vector{T} | "-"
    # Bagnold transport capacity exponent
    e_bagnold::Vector{T} | "-"
end

"Initialize Bagnold transport capacity model parameters"
function TransportCapacityBagnoldParameters(dataset, config, indices)
    lens = lens_input_parameter(
        config,
        "river_water_sediment__bagnold_transport_capacity_coefficient";
        optional = false,
    )
    c_bagnold = ncread(dataset, config, lens; sel = indices, type = Float)
    lens = lens_input_parameter(
        config,
        "river_water_sediment__bagnold_transport_capacity_exponent";
        optional = false,
    )
    e_bagnold = ncread(dataset, config, lens; sel = indices, type = Float)

    tc_parameters =
        TransportCapacityBagnoldParameters(; c_bagnold = c_bagnold, e_bagnold = e_bagnold)

    return tc_parameters
end

"Bagnold river transport capacity model"
@with_kw struct TransportCapacityBagnoldModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityBagnoldParameters{T}
    variables::TransportCapacityModelVariables{T}
end

"Initialize Bagnold river transport capacity model"
function TransportCapacityBagnoldModel(dataset, config, indices)
    n = length(indices)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityBagnoldParameters(dataset, config, indices)
    bc = TransportCapacityBC(n)
    model = TransportCapacityBagnoldModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update Bagnold river transport capacity model for a single timestep"
function update!(model::TransportCapacityBagnoldModel, geometry::RiverGeometry, dt)
    (; q, waterlevel) = model.boundary_conditions
    (; c_bagnold, e_bagnold) = model.parameters
    (; amount) = model.variables

    n = length(q)
    # Note: slope is not used here but this allows for a consistent interface of update! functions
    # Only Bagnold does not use it
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_bagnold(
            q[i],
            waterlevel[i],
            c_bagnold[i],
            e_bagnold[i],
            geometry.width[i],
            geometry.length[i],
            dt,
        )
    end
end

"Engelund and Hansen river transport capacity model parameters"
@with_kw struct TransportCapacityEngelundModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityRiverParameters{T}
    variables::TransportCapacityModelVariables{T}
end

"Initialize Engelund and Hansen river transport capacity model"
function TransportCapacityEngelundModel(dataset, config, indices)
    n = length(indices)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityRiverParameters(dataset, config, indices)
    bc = TransportCapacityBC(n)
    model = TransportCapacityEngelundModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update Engelund and Hansen river transport capacity model for a single timestep"
function update!(model::TransportCapacityEngelundModel, geometry::RiverGeometry, dt)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_engelund(
            q[i],
            waterlevel[i],
            density[i],
            d50[i],
            geometry.width[i],
            geometry.length[i],
            geometry.slope[i],
            dt,
        )
    end
end

"Struct to store Kodatie river transport capacity model parameters"
@get_units @grid_loc @with_kw struct TransportCapacityKodatieParameters{T}
    # Kodatie transport capacity coefficient a
    a_kodatie::Vector{T} | "-"
    # Kodatie transport capacity coefficient b
    b_kodatie::Vector{T} | "-"
    # Kodatie transport capacity coefficient c
    c_kodatie::Vector{T} | "-"
    # Kodatie transport capacity coefficient d
    d_kodatie::Vector{T} | "-"
end

"Initialize Kodatie river transport capacity model parameters"
function TransportCapacityKodatieParameters(dataset, config, indices)
    lens = lens_input_parameter(
        config,
        "river_water_sediment__kodatie_transport_capacity_a-coefficient";
        optional = false,
    )
    a_kodatie = ncread(dataset, config, lens; sel = indices, type = Float)
    lens = lens_input_parameter(
        config,
        "river_water_sediment__kodatie_transport_capacity_b-coefficient";
        optional = false,
    )
    b_kodatie = ncread(dataset, config, lens; sel = indices, type = Float)
    lens = lens_input_parameter(
        config,
        "river_water_sediment__kodatie_transport_capacity_c-coefficient";
        optional = false,
    )
    c_kodatie = ncread(dataset, config, lens; sel = indices, type = Float)
    lens = lens_input_parameter(
        config,
        "river_water_sediment__kodatie_transport_capacity_d-coefficient";
        optional = false,
    )
    d_kodatie = ncread(dataset, config, lens; sel = indices, type = Float)

    tc_parameters = TransportCapacityKodatieParameters(;
        a_kodatie = a_kodatie,
        b_kodatie = b_kodatie,
        c_kodatie = c_kodatie,
        d_kodatie = d_kodatie,
    )

    return tc_parameters
end

"Kodatie river transport capacity model"
@with_kw struct TransportCapacityKodatieModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityKodatieParameters{T}
    variables::TransportCapacityModelVariables{T}
end

"Initialize Kodatie river transport capacity model"
function TransportCapacityKodatieModel(dataset, config, indices)
    n = length(indices)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityKodatieParameters(dataset, config, indices)
    bc = TransportCapacityBC(n)
    model = TransportCapacityKodatieModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update Kodatie river transport capacity model for a single timestep"
function update!(model::TransportCapacityKodatieModel, geometry::RiverGeometry, dt)
    (; q, waterlevel) = model.boundary_conditions
    (; a_kodatie, b_kodatie, c_kodatie, d_kodatie) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_kodatie(
            q[i],
            waterlevel[i],
            a_kodatie[i],
            b_kodatie[i],
            c_kodatie[i],
            d_kodatie[i],
            geometry.width[i],
            geometry.length[i],
            geometry.slope[i],
            dt,
        )
    end
end

"Yang river transport capacity model"
@with_kw struct TransportCapacityYangModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityRiverParameters{T}
    variables::TransportCapacityModelVariables{T}
end

"Initialize Yang river transport capacity model"
function TransportCapacityYangModel(dataset, config, indices)
    n = length(indices)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityRiverParameters(dataset, config, indices)
    bc = TransportCapacityBC(n)
    model = TransportCapacityYangModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update Yang river transport capacity model for a single timestep"
function update!(model::TransportCapacityYangModel, geometry::RiverGeometry, dt)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_yang(
            q[i],
            waterlevel[i],
            density[i],
            d50[i],
            geometry.width[i],
            geometry.length[i],
            geometry.slope[i],
            dt,
        )
    end
end

"Molinas and Wu river transport capacity model"
@with_kw struct TransportCapacityMolinasModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityRiverParameters{T}
    variables::TransportCapacityModelVariables{T}
end

"Initialize Molinas and Wu river transport capacity model"
function TransportCapacityMolinasModel(dataset, config, indices)
    n = length(indices)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityRiverParameters(dataset, config, indices)
    bc = TransportCapacityBC(n)
    model = TransportCapacityMolinasModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

"Update Molinas and Wu river transport capacity model for a single timestep"
function update!(model::TransportCapacityMolinasModel, geometry::RiverGeometry, dt)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_molinas(
            q[i],
            waterlevel[i],
            density[i],
            d50[i],
            geometry.width[i],
            geometry.length[i],
            geometry.slope[i],
            dt,
        )
    end
end