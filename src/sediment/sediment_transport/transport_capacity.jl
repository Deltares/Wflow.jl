abstract type AbstractTransportCapacityModel{T} end

## Total sediment transport capacity structs and functions
@get_units @grid_loc @with_kw struct TransportCapacityModelVariables{T}
    # Total sediment transport capacity
    amount::Vector{T} | "t dt-1"
end

function TransportCapacityModelVariables(n; amount::Vector{T} = fill(mv, n)) where {T}
    return TransportCapacityModelVariables{T}(; amount = amount)
end

@get_units @grid_loc @with_kw struct TransportCapacityBC{T}
    # Discharge
    q::Vector{T} | "m3 s-1"
    # Flow depth
    waterlevel::Vector{T} | "m"
end

function TransportCapacityBC(
    n;
    q::Vector{T} = fill(mv, n),
    waterlevel::Vector{T} = fill(mv, n),
) where {T}
    return TransportCapacityBC{T}(; q = q, waterlevel = waterlevel)
end

function update_boundary_conditions!(
    model::AbstractTransportCapacityModel,
    hydrometeo_forcing::HydrometeoForcing,
    model_type::Symbol,
)
    (; q, waterlevel) = model.boundary_conditions
    (; q_land, waterlevel_land, q_river, waterlevel_river) = hydrometeo_forcing

    if model_type == :land
        @. q = q_land
        @. waterlevel = waterlevel_land
    elseif model_type == :river
        @. q = q_river
        @. waterlevel = waterlevel_river
    end
end

##################### Overland Flow #####################

# Govers parameters for transport capacity models
@get_units @grid_loc @with_kw struct TransportCapacityGoversParameters{T}
    # Drain slope
    slope::Vector{T} | "m m-1"
    # Particle density
    density::Vector{T} | "kg m-3"
    # Govers transport capacity coefficient
    c_govers::Vector{T} | "-"
    # Govers transport capacity exponent
    n_govers::Vector{T} | "-"
end

function TransportCapacityGoversParameters(nc, config, inds)
    slope = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.slope";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    density = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.density";
        sel = inds,
        defaults = 2650.0,
        type = Float,
    )
    c_govers = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.c_govers";
        sel = inds,
        defaults = 0.000505,
        type = Float,
    )
    n_govers = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.n_govers";
        sel = inds,
        defaults = 4.27,
        type = Float,
    )
    tc_parameters = TransportCapacityGoversParameters(;
        slope = slope,
        density = density,
        c_govers = c_govers,
        n_govers = n_govers,
    )

    return tc_parameters
end

@with_kw struct TransportCapacityGoversModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityGoversParameters{T}
    variables::TransportCapacityModelVariables{T}
end

function TransportCapacityGoversModel(nc, config, inds)
    n = length(inds)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityGoversParameters(nc, config, inds)
    bc = TransportCapacityBC(n)
    model = TransportCapacityGoversModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::TransportCapacityGoversModel, width, waterbodies, rivers, ts)
    (; q, waterlevel) = model.boundary_conditions
    (; slope, density, c_govers, n_govers) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_govers(
            q[i],
            waterlevel[i],
            c_govers[i],
            n_govers[i],
            density[i],
            slope[i],
            width[i],
            waterbodies[i],
            rivers[i],
            ts,
        )
    end
end

# Common parameters for transport capacity models
@get_units @grid_loc @with_kw struct TransportCapacityYalinParameters{T}
    # Drain slope
    slope::Vector{T} | "m m-1"
    # Particle density
    density::Vector{T} | "kg m-3"
    # Particle mean diameter
    d50::Vector{T} | "mm"
end

function TransportCapacityYalinParameters(nc, config, inds)
    slope = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.slope";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    density = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.density";
        sel = inds,
        defaults = 2650.0,
        type = Float,
    )
    d50 = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.d50";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    tc_parameters =
        TransportCapacityYalinParameters(; slope = slope, density = density, d50 = d50)

    return tc_parameters
end

@with_kw struct TransportCapacityYalinModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityYalinParameters{T}
    variables::TransportCapacityModelVariables{T}
end

function TransportCapacityYalinModel(nc, config, inds)
    n = length(inds)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityYalinParameters(nc, config, inds)
    bc = TransportCapacityBC(n)
    model = TransportCapacityYalinModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::TransportCapacityYalinModel, width, waterbodies, rivers, ts)
    (; q, waterlevel) = model.boundary_conditions
    (; slope, density, d50) = model.parameters
    (; amount) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_yalin(
            q[i],
            waterlevel[i],
            density[i],
            d50[i],
            slope[i],
            width[i],
            waterbodies[i],
            rivers[i],
            ts,
        )
    end
end

## Total transport capacity with particle differentiation structs and functions
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

# Common parameters for transport capacity models
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

function TransportCapacityYalinDifferentiationParameters(nc, config, inds)
    density = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.density";
        sel = inds,
        defaults = 2650.0,
        type = Float,
    )
    dm_clay = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.dm_clay";
        sel = inds,
        defaults = 2.0,
        type = Float,
    )
    dm_silt = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.dm_silt";
        sel = inds,
        defaults = 10.0,
        type = Float,
    )
    dm_sand = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.dm_sand";
        sel = inds,
        defaults = 200.0,
        type = Float,
    )
    dm_sagg = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.dm_sagg";
        sel = inds,
        defaults = 30.0,
        type = Float,
    )
    dm_lagg = ncread(
        nc,
        config,
        "lateral.land.transport_capacity.parameters.dm_lagg";
        sel = inds,
        defaults = 500.0,
        type = Float,
    )
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

@with_kw struct TransportCapacityYalinDifferentiationModel{T} <:
                AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityYalinDifferentiationParameters{T}
    variables::TransportCapacityYalinDifferentiationModelVariables{T}
end

function TransportCapacityYalinDifferentiationModel(nc, config, inds)
    n = length(inds)
    vars = TransportCapacityYalinDifferentiationModelVariables(n)
    params = TransportCapacityYalinDifferentiationParameters(nc, config, inds)
    bc = TransportCapacityBC(n)
    model = TransportCapacityYalinDifferentiationModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(
    model::TransportCapacityYalinDifferentiationModel,
    geometry::LandParameters,
    waterbodies,
    rivers,
    ts,
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
            ts,
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
            ts,
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
            ts,
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
            ts,
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
            ts,
        )
        amount[i] = clay[i] + silt[i] + sand[i] + sagg[i] + lagg[i]
    end
end

##################### River Flow #####################
@get_units @grid_loc @with_kw struct TransportCapacityRiverParameters{T}
    # Particle density
    density::Vector{T} | "kg m-3"
    # Particle mean diameter
    d50::Vector{T} | "mm"
end

function TransportCapacityRiverParameters(nc, config, inds)
    density = ncread(
        nc,
        config,
        "lateral.river.transport_capacity.parameters.density";
        sel = inds,
        defaults = 2650.0,
        type = Float,
    )
    d50 = ncread(
        nc,
        config,
        "lateral.river.transport_capacity.parameters.d50";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    tc_parameters = TransportCapacityRiverParameters(; density = density, d50 = d50)

    return tc_parameters
end

# Bagnold parameters for transport capacity models
@get_units @grid_loc @with_kw struct TransportCapacityBagnoldParameters{T}
    # Bagnold transport capacity coefficient
    c_bagnold::Vector{T} | "-"
    # Bagnold transport capacity exponent
    e_bagnold::Vector{T} | "-"
end

function TransportCapacityBagnoldParameters(nc, config, inds)
    c_bagnold = ncread(
        nc,
        config,
        "lateral.river.transport_capacity.parameters.c_bagnold";
        sel = inds,
        optional = false,
        type = Float,
    )
    e_bagnold = ncread(
        nc,
        config,
        "lateral.river.transport_capacity.parameters.e_bagnold";
        sel = inds,
        optional = false,
        type = Float,
    )
    tc_parameters =
        TransportCapacityBagnoldParameters(; c_bagnold = c_bagnold, e_bagnold = e_bagnold)

    return tc_parameters
end

@with_kw struct TransportCapacityBagnoldModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityBagnoldParameters{T}
    variables::TransportCapacityModelVariables{T}
end

function TransportCapacityBagnoldModel(nc, config, inds)
    n = length(inds)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityBagnoldParameters(nc, config, inds)
    bc = TransportCapacityBC(n)
    model = TransportCapacityBagnoldModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::TransportCapacityBagnoldModel, geometry::RiverParameters, ts)
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
            ts,
        )
    end
end

# Engelund and Hansen parameters for transport capacity models
@with_kw struct TransportCapacityEngelundModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityRiverParameters{T}
    variables::TransportCapacityModelVariables{T}
end

function TransportCapacityEngelundModel(nc, config, inds)
    n = length(inds)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityRiverParameters(nc, config, inds)
    bc = TransportCapacityBC(n)
    model = TransportCapacityEngelundModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::TransportCapacityEngelundModel, geometry::RiverParameters, ts)
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
            ts,
        )
    end
end

# Kodatie parameters for transport capacity models
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

function TransportCapacityKodatieParameters(nc, config, inds)
    a_kodatie = ncread(
        nc,
        config,
        "lateral.river.transport_capacity.parameters.a_kodatie";
        sel = inds,
        optional = false,
        type = Float,
    )
    b_kodatie = ncread(
        nc,
        config,
        "lateral.river.transport_capacity.parameters.b_kodatie";
        sel = inds,
        optional = false,
        type = Float,
    )
    c_kodatie = ncread(
        nc,
        config,
        "lateral.river.transport_capacity.parameters.c_kodatie";
        sel = inds,
        optional = false,
        type = Float,
    )
    d_kodatie = ncread(
        nc,
        config,
        "lateral.river.transport_capacity.parameters.d_kodatie";
        sel = inds,
        optional = false,
        type = Float,
    )
    tc_parameters = TransportCapacityKodatieParameters(;
        a_kodatie = a_kodatie,
        b_kodatie = b_kodatie,
        c_kodatie = c_kodatie,
        d_kodatie = d_kodatie,
    )

    return tc_parameters
end

@with_kw struct TransportCapacityKodatieModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityKodatieParameters{T}
    variables::TransportCapacityModelVariables{T}
end

function TransportCapacityKodatieModel(nc, config, inds)
    n = length(inds)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityKodatieParameters(nc, config, inds)
    bc = TransportCapacityBC(n)
    model = TransportCapacityKodatieModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::TransportCapacityKodatieModel, geometry::RiverParameters, ts)
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
            ts,
        )
    end
end

# Yang parameters for transport capacity models
@with_kw struct TransportCapacityYangModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityRiverParameters{T}
    variables::TransportCapacityModelVariables{T}
end

function TransportCapacityYangModel(nc, config, inds)
    n = length(inds)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityRiverParameters(nc, config, inds)
    bc = TransportCapacityBC(n)
    model = TransportCapacityYangModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::TransportCapacityYangModel, geometry::RiverParameters, ts)
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
            ts,
        )
    end
end

# Molinas and Wu parameters for transport capacity models
@with_kw struct TransportCapacityMolinasModel{T} <: AbstractTransportCapacityModel{T}
    boundary_conditions::TransportCapacityBC{T}
    parameters::TransportCapacityRiverParameters{T}
    variables::TransportCapacityModelVariables{T}
end

function TransportCapacityMolinasModel(nc, config, inds)
    n = length(inds)
    vars = TransportCapacityModelVariables(n)
    params = TransportCapacityRiverParameters(nc, config, inds)
    bc = TransportCapacityBC(n)
    model = TransportCapacityMolinasModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::TransportCapacityMolinasModel, geometry::RiverParameters, ts)
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
            ts,
        )
    end
end