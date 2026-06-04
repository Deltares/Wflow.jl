abstract type AbstractDemandModel end
abstract type AbstractAllocationModel end

"Land hydrology model with SBM soil model"
@with_kw struct LandHydrologySBM{D<:AbstractDemandModel,A<:AbstractAllocationModel} <:
                AbstractLandModel
    atmospheric_forcing::AtmosphericForcing
    vegetation_parameters::VegetationParameters
    interception::AbstractInterceptionModel
    snow::AbstractSnowModel
    glacier::AbstractGlacierModel
    runoff::AbstractRunoffModel
    soil::SbmSoilModel
    demand::D
    allocation::A
end

"Initialize land hydrology model with SBM soil model"
function LandHydrologySBM(dataset::NCDataset, config::Config, domain::DomainLand)
    (; indices) = domain.network
    dt = Second(config.time.timestepsecs)
    n_cells = length(indices)

    atmospheric_forcing = AtmosphericForcing(; n_cells)
    vegetation_parameters = VegetationParameters(dataset, config, indices)
    if dt >= Hour(23)
        interception =
            GashInterceptionModel(dataset, config, indices, vegetation_parameters)
        @info "Using the Gash interception model since dt >= 23 hours."
    else
        interception = RutterInterceptionModel(vegetation_parameters, n_cells)
        @info "Using the modified Rutter interception model since dt < 23 hours."
    end

    do_snow = config.model.snow__flag
    do_glacier = config.model.glacier__flag
    if do_snow
        snow = SnowHbvModel(dataset, config, indices)
    else
        snow = NoSnowModel(n_cells)
    end
    if do_snow && do_glacier
        glacier_bc = SnowStateBC(; snow_storage=snow.variables.snow_storage)
        glacier = GlacierHbvModel(dataset, config, indices, dt, glacier_bc)
    elseif !do_snow && do_glacier
        @warn string(
            "Glacier processes can be modelled when snow modelling is enabled. To include ",
            "glacier modelling, set `snow__flag` to `true` in the Model section of the TOML file.",
        )
        glacier = NoGlacierModel(n_cells)
    else
        glacier = NoGlacierModel(n_cells)
    end
    runoff = OpenWaterRunoff(; n_cells)

    soil = SbmSoilModel(dataset, config, vegetation_parameters, indices, dt)
    @. vegetation_parameters.rooting_depth =
        min(soil.parameters.soil_thickness * 0.99, vegetation_parameters.rooting_depth)

    if do_water_demand(config)
        allocation = AllocationLandModel(dataset, config, indices)
        demand = DemandModel(dataset, config, indices)
    else
        allocation = NoAllocationLandModel(n_cells)
        demand = NoDemandModel(; n_cells)
    end

    return LandHydrologySBM(;
        atmospheric_forcing,
        vegetation_parameters,
        interception,
        snow,
        glacier,
        runoff,
        soil,
        demand,
        allocation,
    )
end

"Update land hydrology model with SBM soil model for a single timestep"
function update_land_hydrology_model!(
    land_hydrology_model::LandHydrologySBM,
    routing::Routing,
    domain::Domain,
    config::Config,
    dt::Float64,
)
    (; parameters) = domain.land
    (; glacier, snow, interception, runoff, soil, demand, allocation, atmospheric_forcing) =
        land_hydrology_model

    update_interception_model!(interception, atmospheric_forcing, dt)

    update_bc_snow_model!(snow, (; interception))
    update_snow_model!(snow, atmospheric_forcing, dt)

    if config.model.snow_gravitational_transport__flag
        lateral_snow_transport!(snow, domain.land, dt)
    end

    update_glacier_model!(glacier, atmospheric_forcing, dt)

    update_bc_open_water_runoff_model!(
        runoff,
        (; glacier, snow, interception),
        routing,
        domain.river.network,
    )
    update_open_water_runoff_model!(runoff, atmospheric_forcing, parameters, dt)

    if do_water_demand(config)
        (; potential_transpiration) = soil.boundary_conditions
        (; h3_high, h3_low) = soil.parameters
        potential_transpiration .= get_potential_transpiration(interception)
        @. soil.variables.h3 = feddes_h3(h3_high, h3_low, potential_transpiration)
    end
    update_water_demand_model!(demand, soil, dt)
    update_water_allocation_model!(allocation, demand, routing, domain, dt)

    soil_fraction!(soil, glacier, parameters)
    update_bc_soil_model!(
        soil,
        atmospheric_forcing,
        (; interception, runoff, demand, allocation),
        dt,
    )

    update_soil_water_flow!(soil, atmospheric_forcing, (; snow, runoff, demand), config, dt)
    @. soil.variables.actual_evapotranspiration += interception.variables.interception_rate
    return nothing
end

"""
Update the total water storage per cell at the end of a timestep.

# Arguments
- `model`: The land hydrology model `LandHydrologySBM` with the SBM soil model.
- `domain`: Containing shared river and land parameters and network information (e.g. active
    indices).
- `routing`: Containing routing models.
"""
function update_total_water_storage!(
    land_hydrology_model::LandHydrologySBM,
    domain::Domain,
    routing::Routing,
)
    (; overland_flow, river_flow) = routing
    (; interception, snow, glacier, soil, demand) = land_hydrology_model
    (; total_storage, unsaturated_store_depth, saturated_water_depth) = soil.variables

    (; river_fraction, area) = domain.land.parameters
    (; flow_width, flow_length) = domain.river.parameters

    # Set the total storage to zero
    fill!(total_storage, 0)

    # Burn the river routing values
    for (river_cell_idx, index_river) in enumerate(domain.river.network.land_indices)
        total_storage[index_river] = (
            (
                river_flow.variables.h[river_cell_idx] *
                flow_width[river_cell_idx] *
                flow_length[river_cell_idx]
            ) / (area[index_river])
        )
    end

    # Add storage from interception, snow and glacier models
    total_storage .+=
        get_snow_storage(snow) .+ get_snow_water(snow) .+
        get_glacier_store(glacier) .* get_glacier_fraction(glacier) .+
        interception.variables.canopy_storage .+ get_water_depth(demand.paddy)

    # Chunk the data for parallel computing
    n_cells = length(unsaturated_store_depth)
    threaded_foreach(1:n_cells; basesize=1000) do cell_idx
        sub_surface = unsaturated_store_depth[cell_idx] + saturated_water_depth[cell_idx]
        lateral = overland_flow.variables.h[cell_idx] * (1 - river_fraction[cell_idx])

        # Add everything to the total water storage
        total_storage[cell_idx] += sub_surface + lateral
    end
    return nothing
end
