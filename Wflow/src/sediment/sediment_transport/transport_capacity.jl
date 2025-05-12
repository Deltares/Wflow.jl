abstract type AbstractTransportCapacityModel end

"Struct to store total transport capacity model variables"
@with_kw struct TransportCapacityModelVariables
    # Total sediment transport capacity [t dt-1]
    amount::Vector{Float}
end

"Initialize total transport capacity model variables"
function TransportCapacityModelVariables(
    n::Int;
    amount::Vector{Float} = fill(MISSING_VALUE, n),
)
    return TransportCapacityModelVariables(; amount = amount)
end

"Struct to store total transport capacity model boundary conditions"
@with_kw struct TransportCapacityBC
    # Discharge [m³ s⁻¹]
    q::Vector{Float}
    # Flow depth [m]
    waterlevel::Vector{Float}
end

"Initialize total transport capacity model boundary conditions"
function TransportCapacityBC(
    n::Int;
    q::Vector{Float} = fill(MISSING_VALUE, n),
    waterlevel::Vector{Float} = fill(MISSING_VALUE, n),
)
    return TransportCapacityBC(; q = q, waterlevel = waterlevel)
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
@with_kw struct TransportCapacityGoversParameters
    # Particle density [kg m-3]
    density::Vector{Float}
    # Govers transport capacity coefficient [-]
    c_govers::Vector{Float}
    # Govers transport capacity exponent [-]
    n_govers::Vector{Float}
end

"Initialize Govers overland flow transport capacity model parameters"
function TransportCapacityGoversParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "sediment__particle_density")
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
@with_kw struct TransportCapacityGoversModel <: AbstractTransportCapacityModel
    boundary_conditions::TransportCapacityBC
    parameters::TransportCapacityGoversParameters
    variables::TransportCapacityModelVariables
end

"Initialize Govers overland flow transport capacity model"
function TransportCapacityGoversModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
function update!(model::TransportCapacityGoversModel, parameters::LandParameters, dt::Float)
    (; q, waterlevel) = model.boundary_conditions
    (; density, c_govers, n_govers) = model.parameters
    (; amount) = model.variables

    (; slope, flow_width, waterbody_coverage, river_location) = parameters

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_govers(
            q[i],
            waterlevel[i],
            c_govers[i],
            n_govers[i],
            density[i],
            slope[i],
            flow_width[i],
            waterbody_coverage[i],
            river_location[i],
            dt,
        )
    end
end

"Struct to store Yalin overland flow transport capacity model parameters"
@with_kw struct TransportCapacityYalinParameters
    # Particle density [kg m-3]
    density::Vector{Float}
    # Particle mean diameter [mm]
    d50::Vector{Float}
end

"Initialize Yalin overland flow transport capacity model parameters"
function TransportCapacityYalinParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "sediment__particle_density")
    density = ncread(dataset, config, lens; sel = indices, defaults = 2650.0, type = Float)

    lens = lens_input_parameter(config, "land_surface_sediment__median_diameter")
    d50 = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)

    tc_parameters = TransportCapacityYalinParameters(; density = density, d50 = d50)

    return tc_parameters
end

"Yalin overland flow transport capacity model"
@with_kw struct TransportCapacityYalinModel <: AbstractTransportCapacityModel
    boundary_conditions::TransportCapacityBC
    parameters::TransportCapacityYalinParameters
    variables::TransportCapacityModelVariables
end

"Initialize Yalin overland flow transport capacity model"
function TransportCapacityYalinModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
function update!(model::TransportCapacityYalinModel, parameters::LandParameters, dt::Float)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; amount) = model.variables

    (; slope, flow_width, waterbody_coverage, river_location) = parameters

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = transport_capacity_yalin(
            q[i],
            waterlevel[i],
            density[i],
            d50[i],
            slope[i],
            flow_width[i],
            waterbody_coverage[i],
            river_location[i],
            dt,
        )
    end
end

"Struct to store Yalin differentiated overland flow transport capacity model variables"
@with_kw struct TransportCapacityYalinDifferentiationModelVariables
    # Total sediment transport capacity [t dt-1]
    amount::Vector{Float}
    # Transport capacity clay [t dt-1]
    clay::Vector{Float}
    # Transport capacity silt [t dt-1]
    silt::Vector{Float}
    # Transport capacity sand [t dt-1]
    sand::Vector{Float}
    # Transport capacity small aggregates [t dt-1]
    sagg::Vector{Float}
    # Transport capacity large aggregates [t dt-1]
    lagg::Vector{Float}
end

"Initialize Yalin differentiated overland flow transport capacity model variables"
function TransportCapacityYalinDifferentiationModelVariables(
    n::Int;
    amount::Vector{Float} = fill(MISSING_VALUE, n),
    clay::Vector{Float} = fill(MISSING_VALUE, n),
    silt::Vector{Float} = fill(MISSING_VALUE, n),
    sand::Vector{Float} = fill(MISSING_VALUE, n),
    sagg::Vector{Float} = fill(MISSING_VALUE, n),
    lagg::Vector{Float} = fill(MISSING_VALUE, n),
)
    return TransportCapacityYalinDifferentiationModelVariables(;
        amount = amount,
        clay = clay,
        silt = silt,
        sand = sand,
        sagg = sagg,
        lagg = lagg,
    )
end

"Struct to store Yalin differentiated overland flow transport capacity model parameters"
@with_kw struct TransportCapacityYalinDifferentiationParameters
    # Particle density [kg m-3]
    density::Vector{Float}
    # Clay mean diameter [μm]
    dm_clay::Vector{Float}
    # Silt mean diameter [μm]
    dm_silt::Vector{Float}
    # Sand mean diameter [μm]
    dm_sand::Vector{Float}
    # Small aggregates mean diameter [μm]
    dm_sagg::Vector{Float}
    # Large aggregates mean diameter [μm]
    dm_lagg::Vector{Float}
end

"Initialize Yalin differentiated overland flow transport capacity model parameters"
function TransportCapacityYalinDifferentiationParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "sediment__particle_density")
    density = ncread(dataset, config, lens; sel = indices, defaults = 2650.0, type = Float)
    lens = lens_input_parameter(config, "clay__mean_diameter")
    dm_clay = ncread(dataset, config, lens; sel = indices, defaults = 2.0, type = Float)
    lens = lens_input_parameter(config, "silt__mean_diameter")
    dm_silt = ncread(dataset, config, lens; sel = indices, defaults = 10.0, type = Float)
    lens = lens_input_parameter(config, "sand__mean_diameter")
    dm_sand = ncread(dataset, config, lens; sel = indices, defaults = 200.0, type = Float)
    lens = lens_input_parameter(config, "sediment_aggregates~small__mean_diameter")
    dm_sagg = ncread(dataset, config, lens; sel = indices, defaults = 30.0, type = Float)
    lens = lens_input_parameter(config, "sediment_aggregates~large__mean_diameter")
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
@with_kw struct TransportCapacityYalinDifferentiationModel <: AbstractTransportCapacityModel
    boundary_conditions::TransportCapacityBC
    parameters::TransportCapacityYalinDifferentiationParameters
    variables::TransportCapacityYalinDifferentiationModelVariables
end

"Initialize Yalin differentiated overland flow transport capacity model"
function TransportCapacityYalinDifferentiationModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
    parameters::LandParameters,
    dt::Float,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg) = model.parameters
    (; amount, clay, silt, sand, sagg, lagg) = model.variables

    (; slope, flow_width, river_location, waterbody_coverage) = parameters

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
            slope[i],
        )
        clay[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_clay[i],
            slope[i],
            flow_width[i],
            waterbody_coverage[i],
            river_location[i],
            dtot,
            dt,
        )
        silt[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_silt[i],
            slope[i],
            flow_width[i],
            waterbody_coverage[i],
            river_location[i],
            dtot,
            dt,
        )
        sand[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_sand[i],
            slope[i],
            flow_width[i],
            waterbody_coverage[i],
            river_location[i],
            dtot,
            dt,
        )
        sagg[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_sagg[i],
            slope[i],
            flow_width[i],
            waterbody_coverage[i],
            river_location[i],
            dtot,
            dt,
        )
        lagg[i] = transport_capacity_yalin_differentiation(
            q[i],
            waterlevel[i],
            density[i],
            dm_lagg[i],
            slope[i],
            flow_width[i],
            waterbody_coverage[i],
            river_location[i],
            dtot,
            dt,
        )
        amount[i] = clay[i] + silt[i] + sand[i] + sagg[i] + lagg[i]
    end
end

"Struct to store common river transport capacity model parameters"
@with_kw struct TransportCapacityRiverParameters
    # Particle density [kg m-3]
    density::Vector{Float}
    # Particle mean diameter [mm]
    d50::Vector{Float}
end

"Initialize common river transport capacity model parameters"
function TransportCapacityRiverParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "sediment__particle_density")
    density = ncread(dataset, config, lens; sel = indices, defaults = 2650.0, type = Float)
    lens = lens_input_parameter(config, "river_sediment__median_diameter")
    d50 = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)

    tc_parameters = TransportCapacityRiverParameters(; density = density, d50 = d50)

    return tc_parameters
end

"Struct to store Bagnold transport capacity model parameters"
@with_kw struct TransportCapacityBagnoldParameters
    # Bagnold transport capacity coefficient [-]
    c_bagnold::Vector{Float}
    # Bagnold transport capacity exponent [-]
    e_bagnold::Vector{Float}
end

"Initialize Bagnold transport capacity model parameters"
function TransportCapacityBagnoldParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
@with_kw struct TransportCapacityBagnoldModel <: AbstractTransportCapacityModel
    boundary_conditions::TransportCapacityBC
    parameters::TransportCapacityBagnoldParameters
    variables::TransportCapacityModelVariables
end

"Initialize Bagnold river transport capacity model"
function TransportCapacityBagnoldModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
function update!(
    model::TransportCapacityBagnoldModel,
    parameters::RiverParameters,
    dt::Float,
)
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
            parameters.flow_width[i],
            parameters.flow_length[i],
            dt,
        )
    end
end

"Engelund and Hansen river transport capacity model parameters"
@with_kw struct TransportCapacityEngelundModel <: AbstractTransportCapacityModel
    boundary_conditions::TransportCapacityBC
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables
end

"Initialize Engelund and Hansen river transport capacity model"
function TransportCapacityEngelundModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
function update!(
    model::TransportCapacityEngelundModel,
    parameters::RiverParameters,
    dt::Float,
)
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
            parameters.flow_width[i],
            parameters.flow_length[i],
            parameters.slope[i],
            dt,
        )
    end
end

"Struct to store Kodatie river transport capacity model parameters"
@with_kw struct TransportCapacityKodatieParameters
    # Kodatie transport capacity coefficient a [-]
    a_kodatie::Vector{Float}
    # Kodatie transport capacity coefficient b [-]
    b_kodatie::Vector{Float}
    # Kodatie transport capacity coefficient c [-]
    c_kodatie::Vector{Float}
    # Kodatie transport capacity coefficient d [-]
    d_kodatie::Vector{Float}
end

"Initialize Kodatie river transport capacity model parameters"
function TransportCapacityKodatieParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
@with_kw struct TransportCapacityKodatieModel <: AbstractTransportCapacityModel
    boundary_conditions::TransportCapacityBC
    parameters::TransportCapacityKodatieParameters
    variables::TransportCapacityModelVariables
end

"Initialize Kodatie river transport capacity model"
function TransportCapacityKodatieModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
function update!(
    model::TransportCapacityKodatieModel,
    parameters::RiverParameters,
    dt::Float,
)
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
            parameters.flow_width[i],
            parameters.flow_length[i],
            parameters.slope[i],
            dt,
        )
    end
end

"Yang river transport capacity model"
@with_kw struct TransportCapacityYangModel <: AbstractTransportCapacityModel
    boundary_conditions::TransportCapacityBC
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables
end

"Initialize Yang river transport capacity model"
function TransportCapacityYangModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
function update!(model::TransportCapacityYangModel, parameters::RiverParameters, dt::Float)
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
            parameters.flow_width[i],
            parameters.flow_length[i],
            parameters.slope[i],
            dt,
        )
    end
end

"Molinas and Wu river transport capacity model"
@with_kw struct TransportCapacityMolinasModel <: AbstractTransportCapacityModel
    boundary_conditions::TransportCapacityBC
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables
end

"Initialize Molinas and Wu river transport capacity model"
function TransportCapacityMolinasModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
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
function update!(
    model::TransportCapacityMolinasModel,
    parameters::RiverParameters,
    dt::Float,
)
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
            parameters.flow_width[i],
            parameters.flow_length[i],
            parameters.slope[i],
            dt,
        )
    end
end