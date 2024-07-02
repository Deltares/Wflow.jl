@get_units @exchange @grid_type @grid_location @with_kw struct NonIrrigationDemand{T}
    demand_gross::Vector{T}                 # gross water demand [mm Δt⁻¹]
    demand_net::Vector{T}                   # net water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
    returnflow::Vector{T}                   # return flow [mm Δt⁻¹]
end

@get_units @exchange @grid_type @grid_location @with_kw struct NonPaddy{T}
    demand_gross::Vector{T}                     # irrigation gross demand [mm Δt⁻¹] 
    irrigation_efficiency::Vector{T} | "-"      # irrigation efficiency [-]
    maximum_irrigation_depth::Vector{T}         # maximum irrigation depth [mm Δt⁻¹]
    irrigation_areas::Vector{Bool} | "-"        # irrigation areas [-]
    irrigation_trigger::Vector{Bool} | "-"      # irrigation on or off [-]
end

@get_units @exchange @grid_type @grid_location @with_kw struct Paddy{T}
    demand_gross::Vector{T}                     # irrigation gross demand [mm Δt⁻¹] 
    irrigation_efficiency::Vector{T} | "-"      # irrigation efficiency [-]
    maximum_irrigation_depth::Vector{T}         # maximum irrigation depth [mm Δt⁻¹]
    irrigation_areas::Vector{Bool} | "-"        # irrigation areas [-]
    irrigation_trigger::Vector{Bool} | "-"      # irrigation on or off [-]
    h_min::Vector{T} | "mm"                     # minimum required water depth in the irrigated rice field [mm]
    h_opt::Vector{T} | "mm"                     # optimal water depth in the irrigated rice fields [mm]
    h_max::Vector{T} | "mm"                     # water depth when rice field starts spilling water (overflow) [mm]
    h::Vector{T} | "mm"                         # actual water depth in rice field [mm]
end

@get_units @exchange @grid_type @grid_location @with_kw struct WaterAllocationRiver{T}
    act_surfacewater_abst::Vector{T}                    # actual surface water abstraction [mm Δt⁻¹]
    act_surfacewater_abst_vol::Vector{T} | "m3 dt-1"    # actual surface water abstraction [m³ Δt⁻¹]
    available_surfacewater::Vector{T} | "m3"            # available surface water [m³]
    nonirri_returnflow::Vector{T}                       # return flow from non irrigation [mm Δt⁻¹] 
end

@get_units @exchange @grid_type @grid_location @with_kw struct WaterAllocationLand{T}
    irri_demand_gross::Vector{T}                        # irrigation gross demand [mm Δt⁻¹]
    nonirri_demand_gross::Vector{T}                     # non-irrigation gross demand [mm Δt⁻¹]
    total_gross_demand::Vector{T}                       # total gross demand [mm Δt⁻¹]
    frac_sw_used::Vector{T} | "-"                       # fraction surface water used [-]
    areas::Vector{Int} | "-"                            # allocation areas [-]
    surfacewater_demand::Vector{T}                      # demand from surface water [mm Δt⁻¹]
    surfacewater_alloc::Vector{T}                       # allocation from surface water [mm Δt⁻¹]
    act_groundwater_abst::Vector{T}                     # actual groundwater abstraction [mm Δt⁻¹]
    act_groundwater_abst_vol::Vector{T} | "m3 dt-1"     # actual groundwater abstraction [m³ Δt⁻¹]
    available_groundwater::Vector{T} | "m3"             # available groundwater [m³]
    groundwater_demand::Vector{T}                       # demand from groundwater [mm Δt⁻¹]
    groundwater_alloc::Vector{T}                        # allocation from groundwater [mm Δt⁻¹]
    irri_alloc::Vector{T}                               # allocated water for irrigation [mm Δt⁻¹]
    nonirri_alloc::Vector{T}                            # allocated water for non-irrigation [mm Δt⁻¹]
    total_alloc::Vector{T}                              # total allocated water [mm Δt⁻¹]
    nonirri_returnflow::Vector{T}                       # return flow from non irrigation [mm Δt⁻¹]
end

"Return return flow fraction based on gross water demand `demand_gross` and net water demand `demand_net`"
function returnflow_fraction(demand_gross, demand_net)
    fraction = divide(demand_net, demand_gross)
    returnflow_fraction = 1.0 - fraction
    return returnflow_fraction
end

"Initialize water demand for the domestic sector"
function initialize_domestic_demand(nc, config, inds, dt)
    demand_gross =
        ncread(
            nc,
            config,
            "vertical.domestic.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.domestic.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)
    n = length(inds)
    returnflow_f = returnflow_fraction.(demand_gross, demand_net)

    domestic = NonIrrigationDemand{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_f,
        returnflow = fill(Float(0), n),
    )

    return domestic
end

"Initialize water demand for the industry sector"
function initialize_industry_demand(nc, config, inds, dt)
    demand_gross =
        ncread(
            nc,
            config,
            "vertical.industry.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.industry.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)
    n = length(inds)
    returnflow_f = returnflow_fraction.(demand_gross, demand_net)

    industry = NonIrrigationDemand{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_f,
        returnflow = fill(Float(0), n),
    )

    return industry
end

"Initialize water demand for the livestock sector"
function initialize_livestock_demand(nc, config, inds, dt)
    demand_gross =
        ncread(
            nc,
            config,
            "vertical.livestock.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.livestock.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)
    n = length(inds)
    returnflow_f = returnflow_fraction.(demand_gross, demand_net)

    livestock = NonIrrigationDemand{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_f,
        returnflow = fill(Float(0), n),
    )
    return livestock
end

"Initialize paddy (rice) fields for water demand and irrigation computations"
function initialize_paddy(nc, config, inds, dt)
    h_min = ncread(
        nc,
        config,
        "vertical.paddy.h_min";
        sel = inds,
        defaults = 20.0,
        type = Float,
    )
    h_opt = ncread(
        nc,
        config,
        "vertical.paddy.h_opt";
        sel = inds,
        defaults = 50.0,
        type = Float,
    )
    h_max = ncread(
        nc,
        config,
        "vertical.paddy.h_max";
        sel = inds,
        defaults = 80.0,
        type = Float,
    )
    efficiency = ncread(
        nc,
        config,
        "vertical.paddy.irrigation_efficiency";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    areas = ncread(
        nc,
        config,
        "vertical.paddy.irrigation_areas";
        sel = inds,
        optional = false,
        type = Bool,
    )
    irrigation_trigger = ncread(
        nc,
        config,
        "vertical.paddy.irrigation_trigger";
        sel = inds,
        optional = false,
        type = Bool,
    )
    max_irri_depth =
        ncread(
            nc,
            config,
            "vertical.paddy.maximum_irrigation_depth";
            sel = inds,
            defaults = 25.0,
            type = Float,
        ) .* (dt / basetimestep)

    paddy = Paddy{Float}(
        demand_gross = fill(mv, length(inds)),
        irrigation_efficiency = efficiency,
        maximum_irrigation_depth = max_irri_depth,
        irrigation_trigger = irrigation_trigger,
        h_min = h_min,
        h_max = h_max,
        h_opt = h_opt,
        irrigation_areas = areas,
        h = fill(0.0, length(inds)),
    )
    return paddy
end

"Initialize crop (non paddy) fields for water demand and irrigation computations"
function initialize_nonpaddy(nc, config, inds, dt)
    efficiency = ncread(
        nc,
        config,
        "vertical.nonpaddy.irrigation_efficiency";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    areas = ncread(
        nc,
        config,
        "vertical.nonpaddy.irrigation_areas";
        sel = inds,
        defaults = 1,
        optional = false,
        type = Int,
    )
    irrigation_trigger = ncread(
        nc,
        config,
        "vertical.nonpaddy.irrigation_trigger";
        sel = inds,
        defaults = 1,
        optional = false,
        type = Bool,
    )
    max_irri_depth =
        ncread(
            nc,
            config,
            "vertical.nonpaddy.maximum_irrigation_depth";
            sel = inds,
            defaults = 25.0,
            type = Float,
        ) .* (dt / basetimestep)

    nonpaddy = NonPaddy{Float}(
        demand_gross = fill(mv, length(inds)),
        maximum_irrigation_depth = max_irri_depth,
        irrigation_efficiency = efficiency,
        irrigation_areas = areas,
        irrigation_trigger = irrigation_trigger,
    )

    return nonpaddy
end

"Initialize water allocation for the river domain"
function initialize_waterallocation_river(n)
    waterallocation = WaterAllocationRiver(
        act_surfacewater_abst = zeros(Float, n),
        act_surfacewater_abst_vol = zeros(Float, n),
        available_surfacewater = zeros(Float, n),
        nonirri_returnflow = zeros(Float, n),
    )
    return waterallocation
end

"Initialize water allocation for the land domain (`vertical`)"
function initialize_waterallocation_land(nc, config, inds)
    frac_sw_used = ncread(
        nc,
        config,
        "vertical.waterallocation.frac_sw_used";
        sel = inds,
        defaults = 1,
        type = Float,
    )
    areas = ncread(
        nc,
        config,
        "vertical.waterallocation.areas";
        sel = inds,
        defaults = 1,
        type = Int,
    )

    n = length(inds)

    waterallocation = WaterAllocationLand(
        irri_demand_gross = zeros(Float, n),
        nonirri_demand_gross = zeros(Float, n),
        total_gross_demand = zeros(Float, n),
        frac_sw_used = frac_sw_used,
        areas = areas,
        surfacewater_demand = zeros(Float, n),
        surfacewater_alloc = zeros(Float, n),
        act_groundwater_abst = zeros(Float, n),
        act_groundwater_abst_vol = zeros(Float, n),
        available_groundwater = zeros(Float, n),
        groundwater_demand = zeros(Float, n),
        groundwater_alloc = zeros(Float, n),
        irri_alloc = zeros(Float, n),
        nonirri_alloc = zeros(Float, n),
        total_alloc = zeros(Float, n),
        nonirri_returnflow = zeros(Float, n),
    )
    return waterallocation
end

"Return non-irrigation gross demand and update returnflow fraction"
function update_non_irrigation_demand(non_irri::NonIrrigationDemand, i)
    non_irri.returnflow_fraction[i] =
        returnflow_fraction(non_irri.demand_gross[i], non_irri.demand_net[i])
    return non_irri.demand_gross[i]
end

# return zero (gross demand) if non-irrigation sector is not defined
update_non_irrigation_demand(non_irri::Nothing, i) = 0.0

"Update water allocation for river and land domains based on local surface water (river) availability."
function surface_water_allocation_local(land, river, network)
    # maps from the land domain to the internal river domain (linear index), excluding water bodies
    index_river = network.land.index_river_wb
    for i in eachindex(land.waterallocation.surfacewater_demand)
        if index_river[i] > 0.0
            # the available volume is limited by a fixed scaling factor of 0.8 to prevent
            # rivers completely drying out. check for abstraction through inflow (external
            # negative inflow) and adjust available volume.
            if river.inflow[index_river[i]] < 0.0
                inflow = river.inflow[index_river[i]] * land.dt
                available_volume = max(river.volume[index_river[i]] * 0.80 + inflow, 0.0)
            else
                available_volume = river.volume[index_river[i]] * 0.80
            end
            # satisfy surface water demand with available local river volume
            surfacewater_demand_vol =
                land.waterallocation.surfacewater_demand[i] * 0.001 * network.land.area[i]
            abstraction_vol = min(surfacewater_demand_vol, available_volume)
            river.waterallocation.act_surfacewater_abst_vol[index_river[i]] =
                abstraction_vol
            # remaining available surface water and demand 
            river.waterallocation.available_surfacewater[index_river[i]] =
                max(available_volume - abstraction_vol, 0.0)
            abstraction = (abstraction_vol / network.land.area[i]) * 1000.0
            land.waterallocation.surfacewater_demand[i] =
                max(land.waterallocation.surfacewater_demand[i] - abstraction, 0.0)
            # update actual abstraction from river and surface water allocation (land cell)
            river.waterallocation.act_surfacewater_abst[index_river[i]] = abstraction
            land.waterallocation.surfacewater_alloc[i] = abstraction
        end
    end
end

"Update water allocation for river and land domains based on surface water (river) availability for allocation areas."
function surface_water_allocation_area(land, river, network)
    inds_river = network.river.indices_allocation_areas
    inds_land = network.land.indices_allocation_areas
    res_index = network.river.reservoir_index
    lake_index = network.river.lake_index

    m = length(inds_river)
    # loop over allocation areas
    for i = 1:m
        # surface water demand (allocation area)
        sw_demand_vol = 0.0
        for j in inds_land[i]
            sw_demand_vol +=
                land.waterallocation.surfacewater_demand[j] * 0.001 * network.land.area[j]
        end
        # surface water availability (allocation area)
        sw_available = 0.0
        for j in inds_river[i]
            # for reservoir locations use reservoir volume
            if res_index[j] > 0
                k = res_index[j]
                river.waterallocation.available_surfacewater[j] =
                    river.reservoir.volume[k] * 0.98 # limit available reservoir volume
                sw_available += river.waterallocation.available_surfacewater[j]
                # for lake locations use lake volume
            elseif lake_index[j] > 0
                k = lake_index[j]
                river.waterallocation.available_surfacewater[j] =
                    river.lake.storage[k] * 0.98 # limit available lake volume
                sw_available += river.waterallocation.available_surfacewater[j]
                # river volume
            else
                sw_available += river.waterallocation.available_surfacewater[j]
            end
        end
        # total actual surface water abstraction [m3] in an allocation area, minimum of
        # available surface water and demand in an allocation area.
        sw_abstraction = min(sw_available, sw_demand_vol)

        # fraction of available surface water that can be abstracted at allocation area
        # level
        frac_abstract_sw = divide(sw_abstraction, sw_available)
        # fraction of water demand that can be satisfied by available surface water at
        # allocation area level. 
        frac_allocate_sw = divide(sw_abstraction, sw_demand_vol)

        # water abstracted from surface water at each river cell (including reservoir and
        # lake locations).
        for j in inds_river[i]
            river.waterallocation.act_surfacewater_abst_vol[j] +=
                frac_abstract_sw * river.waterallocation.available_surfacewater[j]
            river.waterallocation.act_surfacewater_abst[j] =
                river.waterallocation.act_surfacewater_abst_vol[j] / network.river.area[j]
        end

        # water allocated to each land cell.
        for j in inds_land[i]
            land.waterallocation.surfacewater_alloc[j] +=
                frac_allocate_sw * land.waterallocation.surfacewater_demand[j]
        end
    end
end

"Update water allocation for subsurface domain based on local groundwater availability."
function groundwater_allocation_local(land, groundwater_volume, network)
    for i in eachindex(land.waterallocation.groundwater_demand)
        # land index excluding water bodies
        if network.index_wb[i]
            # groundwater demand based on allocation from surface water.
            land.waterallocation.groundwater_demand[i] = max(
                land.waterallocation.irri_demand_gross[i] +
                land.waterallocation.nonirri_demand_gross[i] -
                land.waterallocation.surfacewater_alloc[i],
                0.0,
            )
            # satisfy groundwater demand with available local groundwater volume
            groundwater_demand_vol =
                land.waterallocation.groundwater_demand[i] * 0.001 * network.area[i]
            available_volume = groundwater_volume[i] * 0.75 # limit available groundwater volume
            abstraction_vol = min(groundwater_demand_vol, available_volume)
            land.waterallocation.act_groundwater_abst_vol[i] = abstraction_vol
            # remaining available groundwater and demand 
            land.waterallocation.available_groundwater[i] =
                max(available_volume - abstraction_vol, 0.0)
            abstraction = (abstraction_vol / network.area[i]) * 1000.0
            land.waterallocation.groundwater_demand[i] =
                max(land.waterallocation.groundwater_demand[i] - abstraction, 0.0)
            # update actual abstraction from groundwater and groundwater allocation (land cell)
            land.waterallocation.act_groundwater_abst[i] = abstraction
            land.waterallocation.groundwater_alloc[i] = abstraction
        end
    end
end

"Update water allocation for subsurface domain based on groundwater availability for allocation areas."
function groundwater_allocation_area(land, network)
    inds_river = network.river.indices_allocation_areas
    inds_land = network.land.indices_allocation_areas
    m = length(inds_river)
    # loop over allocation areas
    for i = 1:m
        # groundwater demand and availability (allocation area)
        gw_demand_vol = 0.0
        gw_available = 0.0
        for j in inds_land[i]
            gw_demand_vol +=
                land.waterallocation.groundwater_demand[j] * 0.001 * network.land.area[j]
            gw_available += land.waterallocation.available_groundwater[j]
        end
        # total actual groundwater abstraction [m3] in an allocation area, minimum of
        # available  groundwater and demand in an allocation area.
        gw_abstraction = min(gw_available, gw_demand_vol)

        # fraction of available groundwater that can be abstracted at allocation area level
        frac_abstract_gw = divide(gw_abstraction, gw_available)
        # fraction of water demand that can be satisfied by available groundwater at
        # allocation area level.
        frac_allocate_gw = divide(gw_abstraction, gw_demand_vol)

        # water abstracted from groundwater and allocated.
        for j in inds_land[i]
            land.waterallocation.act_groundwater_abst_vol[j] +=
                frac_abstract_gw * land.waterallocation.available_groundwater[j]
            land.waterallocation.act_groundwater_abst[j] =
                1000.0 *
                (land.waterallocation.act_groundwater_abst_vol[j] / network.land.area[j])
            land.waterallocation.groundwater_alloc[j] +=
                frac_allocate_gw * land.waterallocation.groundwater_demand[j]
        end
    end
end

"Return and update non-irrigation sector (domestic, livestock, industry) return flow"
function return_flow(non_irri::NonIrrigationDemand, waterallocation)
    for i in eachindex(non_irri.returnflow)
        frac = divide(non_irri.demand_gross[i], waterallocation.nonirri_demand_gross[i])
        allocate = frac * waterallocation.nonirri_alloc[i]
        non_irri.returnflow[i] = non_irri.returnflow_fraction[i] * allocate
    end
    return non_irri.returnflow
end

# return zero (return flow) if non-irrigation sector is not defined
return_flow(non_irri::Nothing, waterallocation) = 0.0

groundwater_volume(gw::LateralSSF) = gw.volume
groundwater_volume(gw) = gw.flow.aquifer.volume

"""
    update_water_allocation(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:Union{SbmModel, SbmGwfModel}}

Update water allocation for model type `sbm` or `sbm_gwf` for a single timestep. First,
surface water abstraction is computed to satisfy local water demand (non-irrigation and
irrigation), and then updated (including lakes and reservoirs) to satisfy the remaining
water demand for allocation areas. Then groundwater abstraction is computed to satisfy the
remaining local water demand, and then updated to satisfy the remaining water demand for
allocation areas. Finally, non-irrigation return flows are updated.
"""
function update_water_allocation(
    model::Model{N,L,V,R,W,T},
) where {N,L,V,R,W,T<:Union{SbmModel,SbmGwfModel}}

    @unpack network, lateral, vertical = model

    river = lateral.river
    index_river = network.land.index_river_wb
    res_index_f = network.river.reservoir_index_f
    lake_index_f = network.river.lake_index_f

    vertical.waterallocation.surfacewater_alloc .= 0.0
    river.waterallocation.act_surfacewater_abst .= 0.0
    river.waterallocation.act_surfacewater_abst_vol .= 0.0
    # total surface water demand for each land cell
    @. vertical.waterallocation.surfacewater_demand =
        vertical.waterallocation.frac_sw_used *
        vertical.waterallocation.nonirri_demand_gross +
        vertical.waterallocation.frac_sw_used * vertical.waterallocation.irri_demand_gross

    # local surface water demand and allocation (river, excluding reservoirs and lakes)
    surface_water_allocation_local(vertical, river, network)
    # surface water demand and allocation for areas
    surface_water_allocation_area(vertical, river, network)

    @. river.abstraction = river.waterallocation.act_surfacewater_abst_vol / vertical.dt

    # for reservoir and lake locations set river abstraction at zero and abstract volume
    # from reservoir and lake, including an update of lake waterlevel
    if !isnothing(river.reservoir)
        @. river.abstraction[res_index_f] = 0.0
        @. river.reservoir.volume -=
            river.waterallocation.act_surfacewater_abst_vol[res_index_f]
    elseif !isnothing(river.lake)
        @. river.abstraction[lake_index_f] = 0.0
        lakes = river.lake
        @. lakes.storage -= river.waterallocation.act_surfacewater_abst_vol[lake_index_f]
        @. lakes.waterlevel =
            waterlevel(lakes.storfunc, lakes.area, lakes.storage, lakes.sh)
    end

    vertical.waterallocation.groundwater_alloc .= 0.0
    vertical.waterallocation.act_groundwater_abst_vol .= 0.0
    vertical.waterallocation.act_groundwater_abst .= 0.0
    # local groundwater demand and allocation
    groundwater_allocation_local(
        vertical,
        groundwater_volume(lateral.subsurface),
        network.land,
    )
    # groundwater demand and allocation for areas
    groundwater_allocation_area(vertical, network)

    # irrigation allocation
    for i in eachindex(vertical.waterallocation.total_alloc)
        vertical.waterallocation.total_alloc[i] =
            vertical.waterallocation.groundwater_alloc[i] +
            vertical.waterallocation.surfacewater_alloc[i]
        frac_irri = divide(
            vertical.waterallocation.irri_demand_gross[i],
            vertical.waterallocation.total_gross_demand[i],
        )
        vertical.waterallocation.irri_alloc[i] =
            frac_irri * vertical.waterallocation.total_alloc[i]
        vertical.waterallocation.nonirri_alloc[i] =
            vertical.waterallocation.total_alloc[i] - vertical.waterallocation.irri_alloc[i]
    end

    # non-irrigation return flows
    returnflow_livestock = return_flow(vertical.livestock, vertical.waterallocation)
    returnflow_domestic = return_flow(vertical.domestic, vertical.waterallocation)
    returnflow_industry = return_flow(vertical.industry, vertical.waterallocation)

    # map non-irrigation return flow to land and river water allocation
    if (
        !isnothing(vertical.livestock) ||
        !isnothing(vertical.domestic) ||
        !isnothing(vertical.industry)
    )
        @. vertical.waterallocation.nonirri_returnflow =
            returnflow_livestock + returnflow_domestic + returnflow_industry

        for i in eachindex(vertical.waterallocation.nonirri_returnflow)
            if index_river[i] > 0.0
                k = index_river[i]
                river.waterallocation.nonirri_returnflow[k] =
                    vertical.waterallocation.nonirri_returnflow[k]
                vertical.waterallocation.nonirri_returnflow[i] = 0.0
            else
                vertical.waterallocation.nonirri_returnflow[i] =
                    vertical.waterallocation.nonirri_returnflow[i]
            end
        end
    end
end
