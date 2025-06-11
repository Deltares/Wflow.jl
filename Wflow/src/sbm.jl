"Land hydrology model with SBM soil model"
@with_kw struct LandHydrologySBM{D, A} <: AbstractLandModel
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
    n = Int(length(indices))

    atmospheric_forcing = AtmosphericForcing(Int(n))
    vegetation_parameters = VegetationParameters(dataset, config, indices)
    if dt >= Hour(23)
        interception =
            GashInterceptionModel(dataset, config, indices, vegetation_parameters)
    else
        interception = RutterInterceptionModel(vegetation_parameters, n)
    end

    do_snow = get(config.model, "snow__flag", false)::Bool
    if do_snow
        snow = SnowHbvModel(dataset, config, indices, dt)
    else
        snow = NoSnowModel()
    end
    do_glacier = get(config.model, "glacier__flag", false)::Bool
    if do_snow && do_glacier
        glacier_bc = SnowStateBC(; snow_storage = snow.variables.snow_storage)
        glacier = GlacierHbvModel(dataset, config, indices, dt, glacier_bc)
    elseif do_snow == false && do_glacier == true
        @warn string(
            "Glacier processes can be modelled when snow modelling is enabled. To include ",
            "glacier modelling, set `snow__flag` to `true` in the Model section of the TOML file.",
        )
        glacier = NoGlacierModel()
    else
        glacier = NoGlacierModel()
    end
    runoff = OpenWaterRunoff(n)

    soil = SbmSoilModel(dataset, config, vegetation_parameters, indices, dt)
    @. vegetation_parameters.rootingdepth =
        min(soil.parameters.soilthickness * 0.99, vegetation_parameters.rootingdepth)

    do_water_demand = haskey(config.model, "water_demand")
    allocation =
        do_water_demand ? AllocationLand(dataset, config, indices) : NoAllocationLand()
    demand = do_water_demand ? Demand(dataset, config, indices, dt) : NoDemand()

    args = (demand, allocation)
    land_hydrology_model = LandHydrologySBM{typeof.(args)...}(;
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
    return land_hydrology_model
end

"Update land hydrology model with SBM soil model for a single timestep"
function update!(
    model::LandHydrologySBM,
    routing::Routing,
    domain::Domain,
    config::Config,
    dt::Float,
)
    do_water_demand = haskey(config.model, "water_demand")::Bool
    (; parameters) = domain.land
    (; glacier, snow, interception, runoff, soil, demand, allocation, atmospheric_forcing) =
        model

    update!(interception, atmospheric_forcing)

    update_boundary_conditions!(snow, (; interception))
    update!(snow, atmospheric_forcing)

    # lateral snow transport
    if get(config.model, "snow_gravitional_transport__flag", false)::Bool
        lateral_snow_transport!(
            snow.variables.snow_storage,
            snow.variables.snow_water,
            parameters.slope,
            domain.land.network,
        )
    end

    update!(glacier, atmospheric_forcing)

    update_boundary_conditions!(
        runoff,
        (; glacier, snow, interception),
        routing,
        domain.river.network,
    )
    update!(runoff, atmospheric_forcing, parameters)

    if do_water_demand
        (; potential_transpiration) = soil.boundary_conditions
        (; h3_high, h3_low) = soil.parameters
        potential_transpiration .= get_potential_transpiration(interception)
        @. soil.variables.h3 = feddes_h3(h3_high, h3_low, potential_transpiration, dt)
    end
    update_water_demand!(demand, soil)
    update_water_allocation!(allocation, demand, routing, domain, dt)

    soil_fraction!(soil, glacier, parameters)
    update_boundary_conditions!(
        soil,
        atmospheric_forcing,
        (; interception, runoff, demand, allocation),
    )

    update!(soil, atmospheric_forcing, (; snow, runoff, demand), config, dt)
    @. soil.variables.actevap += interception.variables.interception_rate
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
    model::LandHydrologySBM,
    domain::Domain,
    routing::Routing,
)
    (; overland_flow, river_flow) = routing
    (; interception, snow, glacier, soil, demand) = model
    (; total_storage, ustoredepth, satwaterdepth) = soil.variables

    (; river_fraction, area) = domain.land.parameters
    (; flow_width, flow_length) = domain.river.parameters

    # Set the total storage to zero
    fill!(total_storage, 0)

    # Burn the river routing values
    for (i, index_river) in enumerate(domain.river.network.land_indices)
        total_storage[index_river] = (
            (river_flow.variables.h_av[i] * flow_width[i] * flow_length[i]) /
            (area[index_river]) * 1000 # Convert to mm
        )
    end

    # Add storage from interception, snow and glacier models
    total_storage .+=
        get_snow_storage(snow) .+ get_snow_water(snow) .+
        get_glacier_store(glacier) .* get_glacier_fraction(glacier) .+
        interception.variables.canopy_storage .+ get_water_depth(demand.paddy)

    # Chunk the data for parallel computing
    AK.foreachindex(total_storage; scheduler = :polyester, min_elems = 1000) do i
        sub_surface = ustoredepth[i] + satwaterdepth[i]
        lateral = (
            overland_flow.variables.h_av[i] * (1 - river_fraction[i]) * 1000 # convert to mm
        )

        # Add everything to the total water storage
        total_storage[i] += (sub_surface + lateral)
    end
    return nothing
end