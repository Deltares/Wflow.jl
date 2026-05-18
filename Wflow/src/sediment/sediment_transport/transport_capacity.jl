abstract type AbstractTransportCapacityModel end

"Struct to store total transport capacity model variables"
@with_kw struct TransportCapacityModelVariables
    n_cells::Int
    # Total sediment transport capacity [t dt-1]
    sediment_transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct to store total transport capacity model boundary conditions"
@with_kw struct TransportCapacityBC
    n_cells::Int
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Flow depth [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Update total transport capacity model boundary conditions"
function update_bc_transport_capacity_model!(
    transport_capacity_model::AbstractTransportCapacityModel,
    hydrological_forcing::HydrologicalForcing,
    model_type::Symbol,
)
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
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
    land_indices_2d::Vector{CartesianIndex{2}},
)
    density = ncread(
        dataset,
        config,
        "sediment__particle_density",
        SoilLossModel;
        sel = land_indices_2d,
    )
    c_govers = ncread(
        dataset,
        config,
        "land_surface_water_sediment__govers_transport_capacity_coefficient",
        SoilLossModel;
        sel = land_indices_2d,
    )
    n_govers = ncread(
        dataset,
        config,
        "land_surface_water_sediment__govers_transport_capacity_exponent",
        SoilLossModel;
        sel = land_indices_2d,
    )
    tc_parameters = TransportCapacityGoversParameters(; density, c_govers, n_govers)

    return tc_parameters
end

"Govers overland flow transport capacity model"
@with_kw struct TransportCapacityGoversModel <: AbstractTransportCapacityModel
    n_cells::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n_cells = n_cells)
    parameters::TransportCapacityGoversParameters
    variables::TransportCapacityModelVariables =
        TransportCapacityModelVariables(; n_cells = n_cells)
end

"Initialize Govers overland flow transport capacity model"
function TransportCapacityGoversModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
)
    n_cells = length(land_indices_2d)
    parameters = TransportCapacityGoversParameters(dataset, config, land_indices_2d)
    transport_capacity_model = TransportCapacityGoversModel(; n_cells, parameters)
    return transport_capacity_model
end

"Update Govers overland flow transport capacity model for a single timestep"
function update_transport_capacity_model!(
    transport_capacity_model::TransportCapacityGoversModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; n_cells) = transport_capacity_model
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
    (; density, c_govers, n_govers) = transport_capacity_model.parameters
    (; sediment_transport_capacity) = transport_capacity_model.variables

    (; slope, flow_width, reservoir_coverage, river_location) = parameters

    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        sediment_transport_capacity[cell_idx] = transport_capacity_govers(
            q[cell_idx],
            waterlevel[cell_idx],
            c_govers[cell_idx],
            n_govers[cell_idx],
            density[cell_idx],
            slope[cell_idx],
            flow_width[cell_idx],
            reservoir_coverage[cell_idx],
            river_location[cell_idx],
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
    land_indices_2d::Vector{CartesianIndex{2}},
)
    density = ncread(
        dataset,
        config,
        "sediment__particle_density",
        SoilLossModel;
        sel = land_indices_2d,
    )
    d50 = ncread(
        dataset,
        config,
        "land_surface_sediment__median_diameter",
        SoilLossModel;
        sel = land_indices_2d,
    )

    tc_parameters = TransportCapacityYalinParameters(; density = density, d50 = d50)

    return tc_parameters
end

"Yalin overland flow transport capacity model"
@with_kw struct TransportCapacityYalinModel <: AbstractTransportCapacityModel
    n_cells::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n_cells = n_cells)
    parameters::TransportCapacityYalinParameters
    variables::TransportCapacityModelVariables =
        TransportCapacityModelVariables(; n_cells = n_cells)
end

"Initialize Yalin overland flow transport capacity model"
function TransportCapacityYalinModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
)
    n_cells = length(land_indices_2d)
    parameters = TransportCapacityYalinParameters(dataset, config, land_indices_2d)
    transport_capacity = TransportCapacityYalinModel(; n_cells, parameters)
    return transport_capacity
end

"Update Yalin overland flow transport capacity model for a single timestep"
function update_transport_capacity_model!(
    transport_capacity_model::TransportCapacityYalinModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; n_cells) = transport_capacity_model
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
    (; density, d50) = transport_capacity_model.parameters
    (; sediment_transport_capacity) = transport_capacity_model.variables

    (; slope, flow_width, reservoir_coverage, river_location) = parameters

    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        sediment_transport_capacity[cell_idx] = transport_capacity_yalin(
            q[cell_idx],
            waterlevel[cell_idx],
            density[cell_idx],
            d50[cell_idx],
            slope[cell_idx],
            flow_width[cell_idx],
            reservoir_coverage[cell_idx],
            river_location[cell_idx],
            dt,
        )
    end
end

"Struct to store Yalin differentiated overland flow transport capacity model variables"
@with_kw struct TransportCapacityYalinDifferentiationModelVariables
    n_cells::Int
    # Total sediment transport capacity [t dt-1]
    sediment_transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Transport capacity clay [t dt-1]
    clay::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Transport capacity silt [t dt-1]
    silt::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Transport capacity sand [t dt-1]
    sand::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Transport capacity small aggregates [t dt-1]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Transport capacity large aggregates [t dt-1]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n_cells)
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
    land_indices_2d::Vector{CartesianIndex{2}},
)
    density = ncread(
        dataset,
        config,
        "sediment__particle_density",
        SoilLossModel;
        sel = land_indices_2d,
    )
    dm_clay =
        ncread(dataset, config, "clay__mean_diameter", SoilLossModel; sel = land_indices_2d)
    dm_silt =
        ncread(dataset, config, "silt__mean_diameter", SoilLossModel; sel = land_indices_2d)
    dm_sand =
        ncread(dataset, config, "sand__mean_diameter", SoilLossModel; sel = land_indices_2d)
    dm_sagg = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter",
        SoilLossModel;
        sel = land_indices_2d,
    )
    dm_lagg = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter",
        SoilLossModel;
        sel = land_indices_2d,
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
    n_cells::Int
    boundary_conditions::TransportCapacityBC = TransportCapacityBC(; n_cells = n_cells)
    parameters::TransportCapacityYalinDifferentiationParameters
    variables::TransportCapacityYalinDifferentiationModelVariables =
        TransportCapacityYalinDifferentiationModelVariables(; n_cells)
end

"Initialize Yalin differentiated overland flow transport capacity model"
function TransportCapacityYalinDifferentiationModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
)
    n_cells = length(land_indices_2d)
    parameters =
        TransportCapacityYalinDifferentiationParameters(dataset, config, land_indices_2d)
    transport_capacity_model =
        TransportCapacityYalinDifferentiationModel(; n_cells, parameters)
    return transport_capacity_model
end

"Update Yalin differentiated overland flow transport capacity model for a single timestep"
function update_transport_capacity_model!(
    transport_capacity_model::TransportCapacityYalinDifferentiationModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; n_cells) = transport_capacity_model
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
    (; density, dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg) =
        transport_capacity_model.parameters
    (; sediment_transport_capacity, clay, silt, sand, sagg, lagg) =
        transport_capacity_model.variables

    (; slope, flow_width, river_location, reservoir_coverage) = parameters

    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        dtot = transportability_yalin_differentiation(
            waterlevel[cell_idx],
            density[cell_idx],
            dm_clay[cell_idx],
            dm_silt[cell_idx],
            dm_sand[cell_idx],
            dm_sagg[cell_idx],
            dm_lagg[cell_idx],
            slope[cell_idx],
        )
        clay[cell_idx] = transport_capacity_yalin_differentiation(
            q[cell_idx],
            waterlevel[cell_idx],
            density[cell_idx],
            dm_clay[cell_idx],
            slope[cell_idx],
            flow_width[cell_idx],
            reservoir_coverage[cell_idx],
            river_location[cell_idx],
            dtot,
            dt,
        )
        silt[cell_idx] = transport_capacity_yalin_differentiation(
            q[cell_idx],
            waterlevel[cell_idx],
            density[cell_idx],
            dm_silt[cell_idx],
            slope[cell_idx],
            flow_width[cell_idx],
            reservoir_coverage[cell_idx],
            river_location[cell_idx],
            dtot,
            dt,
        )
        sand[cell_idx] = transport_capacity_yalin_differentiation(
            q[cell_idx],
            waterlevel[cell_idx],
            density[cell_idx],
            dm_sand[cell_idx],
            slope[cell_idx],
            flow_width[cell_idx],
            reservoir_coverage[cell_idx],
            river_location[cell_idx],
            dtot,
            dt,
        )
        sagg[cell_idx] = transport_capacity_yalin_differentiation(
            q[cell_idx],
            waterlevel[cell_idx],
            density[cell_idx],
            dm_sagg[cell_idx],
            slope[cell_idx],
            flow_width[cell_idx],
            reservoir_coverage[cell_idx],
            river_location[cell_idx],
            dtot,
            dt,
        )
        lagg[cell_idx] = transport_capacity_yalin_differentiation(
            q[cell_idx],
            waterlevel[cell_idx],
            density[cell_idx],
            dm_lagg[cell_idx],
            slope[cell_idx],
            flow_width[cell_idx],
            reservoir_coverage[cell_idx],
            river_location[cell_idx],
            dtot,
            dt,
        )
        sediment_transport_capacity[cell_idx] =
            clay[cell_idx] +
            silt[cell_idx] +
            sand[cell_idx] +
            sagg[cell_idx] +
            lagg[cell_idx]
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
    river_indices_2d::Vector{CartesianIndex{2}},
)
    density = ncread(
        dataset,
        config,
        "sediment__particle_density",
        SoilLossModel;
        sel = river_indices_2d,
    )
    d50 = ncread(
        dataset,
        config,
        "river_sediment__median_diameter",
        SoilLossModel;
        sel = river_indices_2d,
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
    river_indices_2d::Vector{CartesianIndex{2}},
)
    c_bagnold = ncread(
        dataset,
        config,
        "river_water_sediment__bagnold_transport_capacity_coefficient",
        SoilLossModel;
        sel = river_indices_2d,
    )
    e_bagnold = ncread(
        dataset,
        config,
        "river_water_sediment__bagnold_transport_capacity_exponent",
        SoilLossModel;
        sel = river_indices_2d,
    )

    tc_parameters = TransportCapacityBagnoldParameters(; c_bagnold, e_bagnold)

    return tc_parameters
end

"Bagnold river transport capacity model"
@with_kw struct TransportCapacityBagnoldModel <: AbstractTransportCapacityModel
    n_river_cells::Int
    boundary_conditions::TransportCapacityBC =
        TransportCapacityBC(; n_cells = n_river_cells)
    parameters::TransportCapacityBagnoldParameters
    variables::TransportCapacityModelVariables =
        TransportCapacityModelVariables(; n_cells = n_river_cells)
end

"Initialize Bagnold river transport capacity model"
function TransportCapacityBagnoldModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_river_cells = length(river_indices_2d)
    parameters = TransportCapacityBagnoldParameters(dataset, config, river_indices_2d)
    transport_capacity_model = TransportCapacityBagnoldModel(; n_river_cells, parameters)
    return transport_capacity_model
end

"Update Bagnold river transport capacity model for a single timestep"
function update_transport_capacity_model!(
    transport_capacity_model::TransportCapacityBagnoldModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; n_river_cells) = transport_capacity_model
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
    (; c_bagnold, e_bagnold) = transport_capacity_model.parameters
    (; sediment_transport_capacity) = transport_capacity_model.variables

    # Note: slope is not used here but this allows for a consistent interface of update! functions
    # Only Bagnold does not use it
    threaded_foreach(1:n_river_cells; basesize = 1000) do river_cell_idx
        sediment_transport_capacity[river_cell_idx] = transport_capacity_bagnold(
            q[river_cell_idx],
            waterlevel[river_cell_idx],
            c_bagnold[river_cell_idx],
            e_bagnold[river_cell_idx],
            parameters.flow_width[river_cell_idx],
            parameters.flow_length[river_cell_idx],
            dt,
        )
    end
end

"Engelund and Hansen river transport capacity model parameters"
@with_kw struct TransportCapacityEngelundModel <: AbstractTransportCapacityModel
    n_river_cells::Int
    boundary_conditions::TransportCapacityBC =
        TransportCapacityBC(; n_cells = n_river_cells)
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables =
        TransportCapacityModelVariables(; n_cells = n_river_cells)
end

"Initialize Engelund and Hansen river transport capacity model"
function TransportCapacityEngelundModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_river_cells = length(river_indices_2d)
    parameters = TransportCapacityRiverParameters(dataset, config, river_indices_2d)
    transport_capacity_model = TransportCapacityEngelundModel(; n_river_cells, parameters)
    return transport_capacity_model
end

"Update Engelund and Hansen river transport capacity model for a single timestep"
function update_transport_capacity_model!(
    transport_capacity_model::TransportCapacityEngelundModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; n_river_cells) = transport_capacity_model
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
    (; density, d50) = transport_capacity_model.parameters
    (; sediment_transport_capacity) = transport_capacity_model.variables

    threaded_foreach(1:n_river_cells; basesize = 1000) do river_cell_idx
        sediment_transport_capacity[river_cell_idx] = transport_capacity_engelund(
            q[river_cell_idx],
            waterlevel[river_cell_idx],
            density[river_cell_idx],
            d50[river_cell_idx],
            parameters.flow_width[river_cell_idx],
            parameters.flow_length[river_cell_idx],
            parameters.slope[river_cell_idx],
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
    river_indices_2d::Vector{CartesianIndex{2}},
)
    a_kodatie = ncread(
        dataset,
        config,
        "river_water_sediment__kodatie_transport_capacity_a_coefficient",
        SoilLossModel;
        sel = river_indices_2d,
    )
    b_kodatie = ncread(
        dataset,
        config,
        "river_water_sediment__kodatie_transport_capacity_b_coefficient",
        SoilLossModel;
        sel = river_indices_2d,
    )
    c_kodatie = ncread(
        dataset,
        config,
        "river_water_sediment__kodatie_transport_capacity_c_coefficient",
        SoilLossModel;
        sel = river_indices_2d,
    )
    d_kodatie = ncread(
        dataset,
        config,
        "river_water_sediment__kodatie_transport_capacity_d_coefficient",
        SoilLossModel;
        sel = river_indices_2d,
    )

    tc_parameters =
        TransportCapacityKodatieParameters(; a_kodatie, b_kodatie, c_kodatie, d_kodatie)

    return tc_parameters
end

"Kodatie river transport capacity model"
@with_kw struct TransportCapacityKodatieModel <: AbstractTransportCapacityModel
    n_river_cells::Int
    boundary_conditions::TransportCapacityBC =
        TransportCapacityBC(; n_cells = n_river_cells)
    variables::TransportCapacityModelVariables =
        TransportCapacityModelVariables(; n_cells = n_river_cells)
    parameters::TransportCapacityKodatieParameters
end

"Initialize Kodatie river transport capacity model"
function TransportCapacityKodatieModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_river_cells = length(river_indices_2d)
    parameters = TransportCapacityKodatieParameters(dataset, config, river_indices_2d)
    transport_capacity_model = TransportCapacityKodatieModel(; n_river_cells, parameters)
    return transport_capacity_model
end

"Update Kodatie river transport capacity model for a single timestep"
function update_transport_capacity_model!(
    transport_capacity_model::TransportCapacityKodatieModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; n_river_cells) = transport_capacity_model
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
    (; a_kodatie, b_kodatie, c_kodatie, d_kodatie) = transport_capacity_model.parameters
    (; sediment_transport_capacity) = transport_capacity_model.variables

    threaded_foreach(1:n_river_cells; basesize = 1000) do river_cell_idx
        sediment_transport_capacity[river_cell_idx] = transport_capacity_kodatie(
            q[river_cell_idx],
            waterlevel[river_cell_idx],
            a_kodatie[river_cell_idx],
            b_kodatie[river_cell_idx],
            c_kodatie[river_cell_idx],
            d_kodatie[river_cell_idx],
            parameters.flow_width[river_cell_idx],
            parameters.flow_length[river_cell_idx],
            parameters.slope[river_cell_idx],
            dt,
        )
    end
end

"Yang river transport capacity model"
@with_kw struct TransportCapacityYangModel <: AbstractTransportCapacityModel
    n_river_cells::Int
    boundary_conditions::TransportCapacityBC =
        TransportCapacityBC(; n_cells = n_river_cells)
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables =
        TransportCapacityModelVariables(; n_cells = n_river_cells)
end

"Initialize Yang river transport capacity model"
function TransportCapacityYangModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_river_cells = length(river_indices_2d)
    parameters = TransportCapacityRiverParameters(dataset, config, river_indices_2d)
    transport_capacity = TransportCapacityYangModel(; n_river_cells, parameters)
    return transport_capacity
end

"Update Yang river transport capacity model for a single timestep"
function update_transport_capacity_model!(
    transport_capacity_model::TransportCapacityYangModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; n_river_cells) = transport_capacity_model
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
    (; density, d50) = transport_capacity_model.parameters
    (; sediment_transport_capacity) = transport_capacity_model.variables

    threaded_foreach(1:n_river_cells; basesize = 1000) do river_cell_idx
        sediment_transport_capacity[river_cell_idx] = transport_capacity_yang(
            q[river_cell_idx],
            waterlevel[river_cell_idx],
            density[river_cell_idx],
            d50[river_cell_idx],
            parameters.flow_width[river_cell_idx],
            parameters.flow_length[river_cell_idx],
            parameters.slope[river_cell_idx],
            dt,
        )
    end
end

"Molinas and Wu river transport capacity model"
@with_kw struct TransportCapacityMolinasModel <: AbstractTransportCapacityModel
    n_river_cells::Int
    boundary_conditions::TransportCapacityBC =
        TransportCapacityBC(; n_cells = n_river_cells)
    parameters::TransportCapacityRiverParameters
    variables::TransportCapacityModelVariables =
        TransportCapacityModelVariables(; n_cells = n_river_cells)
end

"Initialize Molinas and Wu river transport capacity model"
function TransportCapacityMolinasModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_river_cells = length(river_indices_2d)
    parameters = TransportCapacityRiverParameters(dataset, config, river_indices_2d)
    transport_capacity = TransportCapacityMolinasModel(; n_river_cells, parameters)
    return transport_capacity
end

"Update Molinas and Wu river transport capacity model for a single timestep"
function update_transport_capacity_model!(
    transport_capacity_model::TransportCapacityMolinasModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; n_river_cells) = transport_capacity_model
    (; q, waterlevel) = transport_capacity_model.boundary_conditions
    (; density, d50) = transport_capacity_model.parameters
    (; sediment_transport_capacity) = transport_capacity_model.variables

    threaded_foreach(1:n_river_cells; basesize = 1000) do river_cell_idx
        sediment_transport_capacity[river_cell_idx] = transport_capacity_molinas(
            q[river_cell_idx],
            waterlevel[river_cell_idx],
            density[river_cell_idx],
            d50[river_cell_idx],
            parameters.flow_width[river_cell_idx],
            parameters.flow_length[river_cell_idx],
            parameters.slope[river_cell_idx],
            dt,
        )
    end
end
