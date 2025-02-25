"Land hydrology model with SBM soil model"
@with_kw struct LandHydrologySBM{T, D, A} <: AbstractLandModel
    atmospheric_forcing::AtmosphericForcing{T}
    vegetation_parameter_set::VegetationParameters{T}
    interception::AbstractInterceptionModel{T}
    snow::AbstractSnowModel{T}
    glacier::AbstractGlacierModel{T}
    runoff::AbstractRunoffModel{T}
    soil::SbmSoilModel
    demand::D
    allocation::A
end

"Initialize land hydrology model with SBM soil model"
function LandHydrologySBM(dataset, config, riverfrac, indices)
    dt = Second(config.time.timestepsecs)
    n = length(indices)

    atmospheric_forcing = AtmosphericForcing(n)
    vegetation_parameter_set = VegetationParameters(dataset, config, indices)
    if dt >= Hour(23)
        interception_model =
            GashInterceptionModel(dataset, config, indices, vegetation_parameter_set)
    else
        interception_model = RutterInterceptionModel(vegetation_parameter_set, n)
    end

    modelsnow = get(config.model, "snow", false)::Bool
    if modelsnow
        snow_model = SnowHbvModel(dataset, config, indices, dt)
    else
        snow_model = NoSnowModel{Float64}()
    end
    modelglacier = get(config.model, "glacier", false)::Bool
    if modelsnow && modelglacier
        glacier_bc =
            SnowStateBC{Float64}(; snow_storage = snow_model.variables.snow_storage)
        glacier_model = GlacierHbvModel(dataset, config, indices, dt, glacier_bc)
    elseif modelsnow == false && modelglacier == true
        @warn string(
            "Glacier processes can be modelled when snow modelling is enabled. To include ",
            "glacier modelling, set `snow` to `true` in the Model section of the TOML file.",
        )
        glacier_model = NoGlacierModel{Float64}()
    else
        glacier_model = NoGlacierModel{Float64}()
    end
    runoff_model = OpenWaterRunoff(dataset, config, indices, riverfrac)

    soil_model = SbmSoilModel(dataset, config, vegetation_parameter_set, indices, dt)
    @. vegetation_parameter_set.rootingdepth = min(
        soil_model.parameters.soilthickness * 0.99,
        vegetation_parameter_set.rootingdepth,
    )

    do_water_demand = haskey(config.model, "water_demand")
    allocation =
        do_water_demand ? AllocationLand(dataset, config, indices) :
        NoAllocationLand{Float64}()
    demand = do_water_demand ? Demand(dataset, config, indices, dt) : NoDemand{Float64}()

    args = (demand, allocation)
    land_hydrology_model = LandHydrologySBM{Float64, typeof.(args)...}(;
        atmospheric_forcing = atmospheric_forcing,
        vegetation_parameter_set = vegetation_parameter_set,
        interception = interception_model,
        snow = snow_model,
        glacier = glacier_model,
        runoff = runoff_model,
        soil = soil_model,
        demand = demand,
        allocation = allocation,
    )
    return land_hydrology_model
end

"Update land hydrology model with SBM soil model for a single timestep"
function update!(model::LandHydrologySBM, routing, network, config, dt)
    do_water_demand = haskey(config.model, "water_demand")::Bool
    (; glacier, snow, interception, runoff, soil, demand, allocation, atmospheric_forcing) =
        model

    update!(interception, atmospheric_forcing)

    update_boundary_conditions!(snow, (; interception))
    update!(snow, atmospheric_forcing)

    # lateral snow transport
    if get(config.model, "gravitational_snow_transport", false)::Bool
        lateral_snow_transport!(
            snow.variables.snow_storage,
            snow.variables.snow_water,
            network.land.slope,
            network.land,
        )
    end

    update!(glacier, atmospheric_forcing)

    update_boundary_conditions!(runoff, (; glacier, snow, interception), routing, network)
    update!(runoff, atmospheric_forcing)

    if do_water_demand
        (; potential_transpiration) = soil.boundary_conditions
        (; h3_high, h3_low) = soil.parameters
        potential_transpiration .= get_potential_transpiration(interception)
        @. soil.variables.h3 = feddes_h3(h3_high, h3_low, potential_transpiration, dt)
    end
    update_water_demand!(demand, soil)
    update_water_allocation!(allocation, demand, routing, network, dt)

    soil_fraction!(soil, runoff, glacier)
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

Takes the following parameters:
- model:
    The land hydrology model with the SBM soil model `LandHydrologySBM`
- river_network:
    The indices of the river cells in relation to the active cells, i.e. model.network.index_river
- area:
    Area of the cells acquired from model.network.land.area
- river_routing:
    The river routing struct, i.e. model.routing.river
- land_routing:
    The land routing struct, i.e. model.routing.overland
"""
function update_total_water_storage!(
    model::LandHydrologySBM,
    river_network,
    area,
    river_routing,
    land_routing,
)
    (; interception, snow, glacier, runoff, soil, demand) = model
    (; total_storage, ustoredepth, satwaterdepth) = soil.variables
    (; riverfrac) = runoff.parameters

    # Set the total storage to zero
    fill!(total_storage, 0)

    # Burn the river routing values
    for (i, index_river) in enumerate(river_network)
        total_storage[index_river] = (
            (
                river_routing.variables.h_av[i] *
                river_routing.parameters.flow_width[i] *
                river_routing.parameters.flow_length[i]
            ) / (area[index_river]) * 1000 # Convert to mm
        )
    end

    # Add storage from interception, snow and glacier models
    total_storage .+=
        get_snow_storage(snow) .+ get_snow_water(snow) .+
        get_glacier_store(glacier) .* get_glacier_fraction(glacier) .+
        interception.variables.canopy_storage .+ get_water_depth(demand.paddy)

    # Chunk the data for parallel computing
    n = length(ustoredepth)
    threaded_foreach(1:n; basesize = 1000) do i
        sub_surface = ustoredepth[i] + satwaterdepth[i]
        lateral = (
            land_routing.variables.h_av[i] * (1 - riverfrac[i]) * 1000 # convert to mm
        )

        # Add everything to the total water storage
        total_storage[i] += (sub_surface + lateral)
    end
    return nothing
end