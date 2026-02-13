abstract type AbstractTransportCapacityModel end

"Struct to store total transport capacity model variables"
@with_kw struct TransportCapacityModelVariables
    n::Int
    # Total sediment transport capacity [t dt-1]
    sediment_transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total transport capacity model boundary conditions"
@with_kw struct TransportCapacityBC
    n::Int
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n)
    # Flow depth [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n)
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
    density::Vector{Float64}
    # Govers transport capacity coefficient [-]
    c_govers::Vector{Float64}
    # Govers transport capacity exponent [-]
    n_govers::Vector{Float64}
end

"Initialize Govers overland flow transport capacity model parameters"
function TransportCapacityGoversParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    density = ncread(
        dataset,
        config,
        "sediment__particle_density";
        sel = indices,
        defaults = 2650.0,
        type = Float64,
    )
    c_govers = ncread(
        dataset,
        config,
        "land_surface_water_sediment__govers_transport_capacity_coefficient";
        sel = indices,
        defaults = 0.000505,
        type = Float64,
    )
    n_govers = ncread(
        dataset,
        config,
        "land_surface_water_sediment__govers_transport_capacity_exponent";
        sel = indices,
        defaults = 4.27,
        type = Float64,
    )
    tc_parameters = TransportCapacityGoversParameters(; density, c_govers, n_govers)

    return tc_parameters
end

"Govers overland flow transport capacity model"
@with_kw struct TransportCapacityGoversModel <: AbstractTransportCapacityModel
    n::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n)
    parameters::TransportCapacityGoversParameters
    variables::TransportCapacityModelVariables = TransportCapacityModelVariables(; n)
end

"Initialize Govers overland flow transport capacity model"
function TransportCapacityGoversModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = TransportCapacityGoversParameters(dataset, config, indices)
    model = TransportCapacityGoversModel(; n, parameters)
    return model
end

"Update Govers overland flow transport capacity model for a single timestep"
function update!(
    model::TransportCapacityGoversModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, c_govers, n_govers) = model.parameters
    (; sediment_transport_capacity) = model.variables

    (; slope, flow_width, reservoir_coverage, river_location) = parameters

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        sediment_transport_capacity[i] = transport_capacity_govers(
            q[i],
            waterlevel[i],
            c_govers[i],
            n_govers[i],
            density[i],
            slope[i],
            flow_width[i],
            reservoir_coverage[i],
            river_location[i],
            dt,
        )
    end
end

"Struct to store Yalin overland flow transport capacity model parameters"
@with_kw struct TransportCapacityYalinParameters
    # Particle density [kg m-3]
    density::Vector{Float64}
    # Particle mean diameter [mm]
    d50::Vector{Float64}
end

"Initialize Yalin overland flow transport capacity model parameters"
function TransportCapacityYalinParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    density = ncread(
        dataset,
        config,
        "sediment__particle_density";
        sel = indices,
        defaults = 2650.0,
        type = Float64,
    )
    d50 = ncread(
        dataset,
        config,
        "land_surface_sediment__median_diameter";
        sel = indices,
        defaults = 0.1,
        type = Float64,
    )

    tc_parameters = TransportCapacityYalinParameters(; density = density, d50 = d50)

    return tc_parameters
end

"Yalin overland flow transport capacity model"
@with_kw struct TransportCapacityYalinModel <: AbstractTransportCapacityModel
    n::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n)
    parameters::TransportCapacityYalinParameters
    variables::TransportCapacityModelVariables = TransportCapacityModelVariables(; n)
end

"Initialize Yalin overland flow transport capacity model"
function TransportCapacityYalinModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = TransportCapacityYalinParameters(dataset, config, indices)
    model = TransportCapacityYalinModel(; n, parameters)
    return model
end

"Update Yalin overland flow transport capacity model for a single timestep"
function update!(
    model::TransportCapacityYalinModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; sediment_transport_capacity) = model.variables

    (; slope, flow_width, reservoir_coverage, river_location) = parameters

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        sediment_transport_capacity[i] = transport_capacity_yalin(
            q[i],
            waterlevel[i],
            density[i],
            d50[i],
            slope[i],
            flow_width[i],
            reservoir_coverage[i],
            river_location[i],
            dt,
        )
    end
end

"Struct to store Yalin differentiated overland flow transport capacity model variables"
@with_kw struct TransportCapacityYalinDifferentiationModelVariables
    n::Int
    # Total sediment transport capacity [t dt-1]
    sediment_transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity clay [t dt-1]
    clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity silt [t dt-1]
    silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity sand [t dt-1]
    sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity small aggregates [t dt-1]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity large aggregates [t dt-1]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store Yalin differentiated overland flow transport capacity model parameters"
@with_kw struct TransportCapacityYalinDifferentiationParameters
    # Particle density [kg m-3]
    density::Vector{Float64}
    # Clay mean diameter [μm]
    dm_clay::Vector{Float64}
    # Silt mean diameter [μm]
    dm_silt::Vector{Float64}
    # Sand mean diameter [μm]
    dm_sand::Vector{Float64}
    # Small aggregates mean diameter [μm]
    dm_sagg::Vector{Float64}
    # Large aggregates mean diameter [μm]
    dm_lagg::Vector{Float64}
end

"Initialize Yalin differentiated overland flow transport capacity model parameters"
function TransportCapacityYalinDifferentiationParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    density = ncread(
        dataset,
        config,
        "sediment__particle_density";
        sel = indices,
        defaults = 2650.0,
        type = Float64,
    )
    dm_clay = ncread(
        dataset,
        config,
        "clay__mean_diameter";
        sel = indices,
        defaults = 2.0,
        type = Float64,
    )
    dm_silt = ncread(
        dataset,
        config,
        "silt__mean_diameter";
        sel = indices,
        defaults = 10.0,
        type = Float64,
    )
    dm_sand = ncread(
        dataset,
        config,
        "sand__mean_diameter";
        sel = indices,
        defaults = 200.0,
        type = Float64,
    )
    dm_sagg = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter";
        sel = indices,
        defaults = 30.0,
        type = Float64,
    )
    dm_lagg = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter";
        sel = indices,
        defaults = 500.0,
        type = Float64,
    )

    tc_parameters = TransportCapacityYalinDifferentiationParameters(;
        density,
        dm_clay,
        dm_silt,
        dm_sand,
        dm_sagg,
        dm_lagg,
    )

    return tc_parameters
end

"Yalin differentiated overland flow transport capacity model"
@with_kw struct TransportCapacityYalinDifferentiationModel <: AbstractTransportCapacityModel
    n::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n)
    parameters::TransportCapacityYalinDifferentiationParameters
    variables::TransportCapacityYalinDifferentiationModelVariables =
        TransportCapacityYalinDifferentiationModelVariables(; n)
end

"Initialize Yalin differentiated overland flow transport capacity model"
function TransportCapacityYalinDifferentiationModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = TransportCapacityYalinDifferentiationParameters(dataset, config, indices)
    model = TransportCapacityYalinDifferentiationModel(; n, parameters)
    return model
end

"Update Yalin differentiated overland flow transport capacity model for a single timestep"
function update!(
    model::TransportCapacityYalinDifferentiationModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg) = model.parameters
    (; sediment_transport_capacity, clay, silt, sand, sagg, lagg) = model.variables

    (; slope, flow_width, river_location, reservoir_coverage) = parameters

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
            reservoir_coverage[i],
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
            reservoir_coverage[i],
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
            reservoir_coverage[i],
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
            reservoir_coverage[i],
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
            reservoir_coverage[i],
            river_location[i],
            dtot,
            dt,
        )
        sediment_transport_capacity[i] = clay[i] + silt[i] + sand[i] + sagg[i] + lagg[i]
    end
end

"Struct to store common river transport capacity model parameters"
@with_kw struct TransportCapacityRiverParameters
    # Particle density [kg m-3]
    density::Vector{Float64}
    # Particle mean diameter [mm]
    d50::Vector{Float64}
end

"Initialize common river transport capacity model parameters"
function TransportCapacityRiverParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    density = ncread(
        dataset,
        config,
        "sediment__particle_density";
        sel = indices,
        defaults = 2650.0,
        type = Float64,
    )
    d50 = ncread(
        dataset,
        config,
        "river_sediment__median_diameter";
        sel = indices,
        defaults = 0.1,
        type = Float64,
    )

    tc_parameters = TransportCapacityRiverParameters(; density, d50)

    return tc_parameters
end

"Struct to store Bagnold transport capacity model parameters"
@with_kw struct TransportCapacityBagnoldParameters
    # Bagnold transport capacity coefficient [-]
    c_bagnold::Vector{Float64}
    # Bagnold transport capacity exponent [-]
    e_bagnold::Vector{Float64}
end

"Initialize Bagnold transport capacity model parameters"
function TransportCapacityBagnoldParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    c_bagnold = ncread(
        dataset,
        config,
        "river_water_sediment__bagnold_transport_capacity_coefficient";
        optional = false,
        sel = indices,
        type = Float64,
    )
    e_bagnold = ncread(
        dataset,
        config,
        "river_water_sediment__bagnold_transport_capacity_exponent";
        optional = false,
        sel = indices,
        type = Float64,
    )

    tc_parameters = TransportCapacityBagnoldParameters(; c_bagnold, e_bagnold)

    return tc_parameters
end

"Bagnold river transport capacity model"
@with_kw struct TransportCapacityBagnoldModel <: AbstractTransportCapacityModel
    n::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n)
    parameters::TransportCapacityBagnoldParameters
    variables::TransportCapacityModelVariables = TransportCapacityModelVariables(; n)
end

"Initialize Bagnold river transport capacity model"
function TransportCapacityBagnoldModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = TransportCapacityBagnoldParameters(dataset, config, indices)
    model = TransportCapacityBagnoldModel(; n, parameters)
    return model
end

"Update Bagnold river transport capacity model for a single timestep"
function update!(
    model::TransportCapacityBagnoldModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; q, waterlevel) = model.boundary_conditions
    (; c_bagnold, e_bagnold) = model.parameters
    (; sediment_transport_capacity) = model.variables

    n = length(q)
    # Note: slope is not used here but this allows for a consistent interface of update! functions
    # Only Bagnold does not use it
    threaded_foreach(1:n; basesize = 1000) do i
        sediment_transport_capacity[i] = transport_capacity_bagnold(
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
    n::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n)
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables = TransportCapacityModelVariables(; n)
end

"Initialize Engelund and Hansen river transport capacity model"
function TransportCapacityEngelundModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = TransportCapacityRiverParameters(dataset, config, indices)
    model = TransportCapacityEngelundModel(; n, parameters)
    return model
end

"Update Engelund and Hansen river transport capacity model for a single timestep"
function update!(
    model::TransportCapacityEngelundModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; sediment_transport_capacity) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        sediment_transport_capacity[i] = transport_capacity_engelund(
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
    a_kodatie::Vector{Float64}
    # Kodatie transport capacity coefficient b [-]
    b_kodatie::Vector{Float64}
    # Kodatie transport capacity coefficient c [-]
    c_kodatie::Vector{Float64}
    # Kodatie transport capacity coefficient d [-]
    d_kodatie::Vector{Float64}
end

"Initialize Kodatie river transport capacity model parameters"
function TransportCapacityKodatieParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    a_kodatie = ncread(
        dataset,
        config,
        "river_water_sediment__kodatie_transport_capacity_a_coefficient";
        optional = false,
        sel = indices,
        type = Float64,
    )
    b_kodatie = ncread(
        dataset,
        config,
        "river_water_sediment__kodatie_transport_capacity_b_coefficient";
        optional = false,
        sel = indices,
        type = Float64,
    )
    c_kodatie = ncread(
        dataset,
        config,
        "river_water_sediment__kodatie_transport_capacity_c_coefficient";
        optional = false,
        sel = indices,
        type = Float64,
    )
    d_kodatie = ncread(
        dataset,
        config,
        "river_water_sediment__kodatie_transport_capacity_d_coefficient";
        optional = false,
        sel = indices,
        type = Float64,
    )

    tc_parameters =
        TransportCapacityKodatieParameters(; a_kodatie, b_kodatie, c_kodatie, d_kodatie)

    return tc_parameters
end

"Kodatie river transport capacity model"
@with_kw struct TransportCapacityKodatieModel <: AbstractTransportCapacityModel
    n::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n)
    variables::TransportCapacityModelVariables = TransportCapacityModelVariables(; n)
    parameters::TransportCapacityKodatieParameters
end

"Initialize Kodatie river transport capacity model"
function TransportCapacityKodatieModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = TransportCapacityKodatieParameters(dataset, config, indices)
    model = TransportCapacityKodatieModel(; n, parameters)
    return model
end

"Update Kodatie river transport capacity model for a single timestep"
function update!(
    model::TransportCapacityKodatieModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; q, waterlevel) = model.boundary_conditions
    (; a_kodatie, b_kodatie, c_kodatie, d_kodatie) = model.parameters
    (; sediment_transport_capacity) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        sediment_transport_capacity[i] = transport_capacity_kodatie(
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
    n::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n)
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables = TransportCapacityModelVariables(; n)
end

"Initialize Yang river transport capacity model"
function TransportCapacityYangModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = TransportCapacityRiverParameters(dataset, config, indices)
    model = TransportCapacityYangModel(; n, parameters)
    return model
end

"Update Yang river transport capacity model for a single timestep"
function update!(
    model::TransportCapacityYangModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; sediment_transport_capacity) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        sediment_transport_capacity[i] = transport_capacity_yang(
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
    n::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n)
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables = TransportCapacityModelVariables(; n)
end

"Initialize Molinas and Wu river transport capacity model"
function TransportCapacityMolinasModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = TransportCapacityRiverParameters(dataset, config, indices)
    model = TransportCapacityMolinasModel(; n, parameters)
    return model
end

"Update Molinas and Wu river transport capacity model for a single timestep"
function update!(
    model::TransportCapacityMolinasModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; q, waterlevel) = model.boundary_conditions
    (; density, d50) = model.parameters
    (; sediment_transport_capacity) = model.variables

    n = length(q)
    threaded_foreach(1:n; basesize = 1000) do i
        sediment_transport_capacity[i] = transport_capacity_molinas(
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
