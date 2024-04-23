@get_units @exchange @grid_type @grid_location @with_kw struct Industry{T}
    demand_gross::Vector{T}                 # gross industry water demand [mm Δt⁻¹]
    demand_net::Vector{T}                   # net industry water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
    returnflow::Vector{T}                   # return flow [mm Δt⁻¹]
end

@get_units @exchange @grid_type @grid_location @with_kw struct Domestic{T}
    demand_gross::Vector{T}                 # gross domestic water demand [mm Δt⁻¹]
    demand_net::Vector{T}                   # net domestic water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
    returnflow::Vector{T}                   # return flow [mm Δt⁻¹]
end

@get_units @exchange @grid_type @grid_location @with_kw struct Livestock{T}
    demand_gross::Vector{T}                 # gross livestock water demand [mm Δt⁻¹]
    demand_net::Vector{T}                   # net livestock water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
    returnflow::Vector{T}                   # return flow [mm Δt⁻¹]
end

@get_units @exchange @grid_type @grid_location @with_kw struct NonPaddy{T}
    demand_gross::Vector{T}                 # irrigation gross demand [mm Δt⁻¹] 
    irrigation_efficiency::Vector{T} | "-"  # irrigation efficiency [-]
    irrigation_areas::Vector{Bool} | "-"    # irrigation areas [-]
    irrigation_trigger::Vector{Bool} | "-"  # irrigation on or off [-]
end

@get_units @exchange @grid_type @grid_location @with_kw struct Paddy{T}
    demand_gross::Vector{T}                 # irrigation gross demand [mm Δt⁻¹] 
    irrigation_efficiency::Vector{T} | "-"  # irrigation efficiency [-]
    irrigation_areas::Vector{Bool} | "-"    # irrigation areas [-]
    irrigation_trigger::Vector{Bool} | "-"  # irrigation on or off [-]
    h_min::Vector{T} | "mm"                 # minimum required water depth in the irrigated rice field [mm]
    h_opt::Vector{T} | "mm"                 # optimal water depth in the irrigated rice fields [mm]
    h_max::Vector{T} | "mm"                 # water depth when rice field starts spilling water (overflow) [mm]
    h::Vector{T} | "mm"                     # actual water depth in rice field [mm]
end

@get_units @exchange @grid_type @grid_location @with_kw struct WaterAllocationRiver{T}
    act_surfacewater_abst::Vector{T}                    # actual surface water abstraction [mm Δt⁻¹]
    act_surfacewater_abst_vol::Vector{T} | "m3 dt-1"    # actual surface water abstraction [m³ Δt⁻¹]
    available_surfacewater::Vector{T} | "m3"            # available surface water [m³]
    nonirri_returnflow::Vector{T}                       # return flow from non irrigation [mm Δt⁻¹] 
end

@get_units @exchange @grid_type @grid_location @with_kw struct WaterAllocationLand{T}
    irri_demand_gross::Vector{T}                        # irrigation gross demand [mm Δt⁻¹]
    nonirri_demand_gross::Vector{T}                     # non irrigation (industry and domestic) gross demand [mm Δt⁻¹]
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

"Set return flow fraction based on gross water demand `demand_gross` and net water demand `demand_net`"
function set_returnflow_fraction(returnflow_fraction, demand_gross, demand_net)
    for i in eachindex(demand_gross)
        fraction = divide(demand_net[i], demand_gross[i])
        returnflow_fraction[i] = 1.0 - fraction
    end
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
    returnflow_fraction = set_returnflow_fraction(fill(mv, n), demand_gross, demand_net)

    domestic = Domestic{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
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
    returnflow_fraction = set_returnflow_fraction(fill(mv, n), demand_gross, demand_net)

    industry = Industry{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
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
    returnflow_fraction = set_returnflow_fraction(fill(mv, n), demand_gross, demand_net)

    livestock = Livestock{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
        returnflow = fill(Float(0), n),
    )
    return livestock
end

"Initialize paddy (rice) fields for water demand and irrigation computations"
function initialize_paddy(nc, config, inds)
    h_min =
        ncread(nc, config, "vertical.paddy.h_min"; sel = inds, defaults = 0.0, type = Float)
    h_opt =
        ncread(nc, config, "vertical.paddy.h_opt"; sel = inds, defaults = 0.0, type = Float)
    h_max =
        ncread(nc, config, "vertical.paddy.h_max"; sel = inds, defaults = 0.0, type = Float)
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
        type = Int,
    )
    irrigation_trigger = ncread(
        nc,
        config,
        "vertical.paddy.irrigation_trigger";
        sel = inds,
        optional = false,
        type = Bool,
    )

    paddy = Paddy{Float}(
        demand_gross = fill(mv, length(inds)),
        irrigation_efficiency = efficiency,
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
function initialize_nonpaddy(nc, config, inds)
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

    nonpaddy = NonPaddy{Float}(
        demand_gross = fill(mv, length(inds)),
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
        #optional = false,
        type = Float,
        #fill = 0,
    )
    areas = ncread(
        nc,
        config,
        "vertical.waterallocation.areas";
        sel = inds,
        defaults = 1,
        #optional = false,
        type = Int,
        #fill = 0,
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

"""
    update_water_demand(sbm::SBM)

Update water demand for vertical `SBM` concept for a single timestep. Water demand is
computed for sectors `industry`, `domestic` and `livestock`, and `paddy` rice fields and
`nonpaddy` (other crop) fields. 

Gross water demand for irrigation `irri_demand_gross` and non-irrigation
`nonirri_demand_gross`, and total gross water demand `total_gross_demand` are updated as
part of `SBM` water allocation (`waterallocation`).
"""
function update_water_demand(sbm::SBM)
    for i = 1:sbm.n
        industry_dem = isnothing(sbm.industry) ? 0.0 : sbm.industry.demand_gross[i]
        domestic_dem = isnothing(sbm.domestic) ? 0.0 : sbm.domestic.demand_gross[i]
        livestock_dem = isnothing(sbm.livestock) ? 0.0 : sbm.livestock.demand_gross[i]

        irri_dem_gross = 0.0
        if !isnothing(sbm.nonpaddy) && sbm.nonpaddy.irrigation_areas[i]
            if sbm.nonpaddy.irrigation_trigger[i]
                usl, _ = set_layerthickness(sbm.zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
                for k = 1:sbm.n_unsatlayers[i]
                    rootfrac = min(
                        1.0,
                        (max(0.0, sbm.rootingdepth[i] - sbm.sumlayers[i][k]) / usl[k]),
                    )
                    vwc_fc = vwc_brooks_corey(
                        -100.0,
                        sbm.hb[i],
                        sbm.θₛ[i],
                        sbm.θᵣ[i],
                        sbm.c[i][k],
                    )
                    vwc_h3 = vwc_brooks_corey(
                        sbm.h3[i],
                        sbm.hb[i],
                        sbm.θₛ[i],
                        sbm.θᵣ[i],
                        sbm.c[i][k],
                    )
                    depletion = (vwc_fc * usl[k]) - sbm.ustorelayerdepth[i][k]
                    depletion *= rootfrac
                    raw = (vwc_fc - vwc_h3) * usl[k] # readily available water
                    raw *= rootfrac
                    if depletion >= raw
                        irri_dem_gross += depletion
                    end
                end
                # limit irrigation demand to infiltration capacity    
                infiltration_capacity =
                    sbm.soilinfredu[i] * (sbm.infiltcappath[i] + sbm.infiltcapsoil[i])
                irri_dem_gross = min(irri_dem_gross, infiltration_capacity)
                irri_dem_gross /= sbm.nonpaddy.irrigation_efficiency[i]
            else
                irri_dem_gross = 0.0
            end
        elseif !isnothing(sbm.paddy) && sbm.paddy.irrigation_areas[i]
            if sbm.paddy.irrigation_trigger[i]
                irr_depth_paddy =
                    sbm.paddy.h[i] < sbm.paddy.h_min[i] ?
                    (sbm.paddy.h_opt[i] - sbm.paddy.h[i]) : 0.0
                irri_dem_gross += irr_depth_paddy / sbm.paddy.irrigation_efficiency[i]
            end
        end
        sbm.waterallocation.irri_demand_gross[i] = irri_dem_gross
        sbm.waterallocation.nonirri_demand_gross[i] =
            industry_dem + domestic_dem + livestock_dem
        sbm.waterallocation.total_gross_demand[i] =
            irri_dem_gross + industry_dem + domestic_dem + livestock_dem
    end
end

"Update water allocation for river and land domains based on local surface water (river) availability."
function surface_water_allocation_local(land, river, network)
    index_river = network.land.index_river_wb
    for i in eachindex(land.waterallocation.surfacewater_demand)
        if index_river[i] > 0.0
            # check for abstraction through inflow (external negative inflow and adjust
            # available volume
            if river.inflow[index_river[i]] < 0.0
                inflow = river.inflow[index_river[i]] * land.dt
                available_volume = max(river.volume[index_river[i]] * 0.80 + inflow, 0.0)
            else
                available_volume = river.volume[index_river[i]] * 0.80
            end
            surfacewater_demand_vol =
                land.waterallocation.surfacewater_demand[i] * 0.001 * network.land.area[i]
            abstraction_vol = min(surfacewater_demand_vol, available_volume)
            river.waterallocation.act_surfacewater_abst_vol[index_river[i]] =
                abstraction_vol
            river.waterallocation.available_surfacewater[index_river[i]] =
                available_volume - abstraction_vol
            abstraction = (abstraction_vol / network.land.area[i]) * 1000.0
            land.waterallocation.surfacewater_demand[i] -= abstraction
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
            if res_index[j] > 0
                k = res_index[j]
                river.waterallocation.available_surfacewater[j] =
                    river.reservoir.volume[k] * 0.98
                sw_available += river.waterallocation.available_surfacewater[j]
            elseif lake_index[j] > 0
                k = lake_index[j]
                river.waterallocation.available_surfacewater[j] =
                    river.lake.storage[k] * 0.98
                sw_available += river.waterallocation.available_surfacewater[j]
            else
                sw_available += river.waterallocation.available_surfacewater[j]
            end
        end
        sw_abstraction = min(sw_available, sw_demand_vol)

        frac_abstract_sw = divide(sw_abstraction, sw_available)
        frac_allocate_sw = divide(sw_abstraction, sw_demand_vol)

        for j in inds_river[i]
            river.waterallocation.act_surfacewater_abst_vol[j] +=
                frac_abstract_sw * river.waterallocation.available_surfacewater[j]
            river.waterallocation.act_surfacewater_abst[j] =
                river.waterallocation.act_surfacewater_abst_vol[j] / network.river.area[j]
        end

        for j in inds_land[i]
            land.waterallocation.surfacewater_alloc[j] +=
                frac_allocate_sw * land.waterallocation.surfacewater_demand[j]
        end
    end
end

"Update water allocation for subsurface domain based on local groundwater availability."
function groundwater_allocation_local(land, groundwater, area)
    for i in eachindex(land.waterallocation.groundwater_demand)
        land.waterallocation.groundwater_demand[i] =
            land.waterallocation.irri_demand_gross[i] +
            land.waterallocation.nonirri_demand_gross[i] -
            land.waterallocation.surfacewater_alloc[i]
        groundwater_demand_vol =
            land.waterallocation.groundwater_demand[i] * 0.001 * area[i]
        available_volume = groundwater.volume[i] * 0.75
        abstraction_vol = min(groundwater_demand_vol, available_volume)
        land.waterallocation.available_groundwater[i] = available_volume - abstraction_vol
        abstraction = (abstraction_vol / area[i]) * 1000.0
        land.waterallocation.groundwater_demand[i] -= abstraction
        land.waterallocation.act_groundwater_abst[i] = abstraction
        land.waterallocation.groundwater_alloc[i] = abstraction
    end
end

"Update water allocation for subsurface domain based on groundwater availability for allocation areas."
function groundwater_allocation_area(land, network)
    inds_river = network.river.indices_allocation_areas
    inds_land = network.land.indices_allocation_areas
    m = length(inds_river)
    for i = 1:m
        gw_demand_vol = 0.0
        gw_available = 0.0
        for j in inds_land[i]
            gw_demand_vol +=
                land.waterallocation.groundwater_demand[j] * 0.001 * network.land.area[j]
            gw_available += land.waterallocation.available_groundwater[j]
        end
        gw_abstraction = min(gw_available, gw_demand_vol)

        frac_abstract_gw = divide(gw_abstraction, gw_available)
        frac_allocate_gw = divide(gw_abstraction, gw_demand_vol)

        for j in inds_land[i]
            land.waterallocation.act_groundwater_abst_vol[j] +=
                frac_abstract_gw * land.waterallocation.available_groundwater[j]
            land.waterallocation.act_groundwater_abst[j] =
                land.waterallocation.act_groundwater_abst_vol[j] / network.land.area[j]
            land.waterallocation.groundwater_alloc[j] +=
                frac_allocate_gw * land.waterallocation.groundwater_demand[j]
        end
    end
end

"Return livestock sector `returnflow`"
function livestock_return_flow(land)
    if isnothing(land.livestock)
        return 0.0
    else
        for i in eachindex(land.livestock.returnflow)
            frac_livestock = divide(
                land.livestock.demand_gross[i],
                land.waterallocation.nonirri_demand_gross[i],
            )
            livestock_alloc = frac_livestock * land.waterallocation.nonirri_alloc[i]
            land.livestock.returnflow[i] =
                land.livestock.returnflow_fraction[i] * livestock_alloc
        end
        return land.livestock.returnflow
    end
end

"Return domestic sector `returnflow`"
function domestic_return_flow(land)
    if isnothing(land.domestic)
        return 0.0
    else
        for i in eachindex(land.domestic.returnflow)
            frac_domestic = divide(
                land.domestic.demand_gross[i],
                land.waterallocation.nonirri_demand_gross[i],
            )
            domestic_alloc = frac_domestic * land.waterallocation.nonirri_alloc[i]
            land.domestic.returnflow[i] =
                land.domestic.returnflow_fraction[i] * domestic_alloc
        end
        return land.domestic.returnflow
    end
end

"Return industry sector `returnflow`"
function industry_return_flow(land)
    if isnothing(land.industry)
        return 0.0
    else
        for i in eachindex(land.industry.returnflow)
            frac_industry = divide(
                land.industry.demand_gross[i],
                land.waterallocation.nonirri_demand_gross[i],
            )
            industry_alloc = frac_industry * land.waterallocation.nonirri_alloc[i]
            land.industry.returnflow[i] =
                land.industry.returnflow_fraction[i] * industry_alloc
        end
        return land.industry.returnflow
    end
end

"""
    update_water_allocation(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

Update water allocation for model type `sbm` for a single timestep. First, surface water
abstraction is computed to satisfy local water demand (non-irrigation and irrigation), and
then updated (including lakes and reservoirs) to satisfy the remaining water demand for
allocation areas. Then groundwater abstraction is computed to satisfy the remaining local
water demand, and then updated to satisfy the remaining water demand for allocation areas.
Finally, non-irrigation return flows are updated.
"""
function update_water_allocation(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

    @unpack network, lateral, vertical = model

    river = lateral.river
    index_river = network.land.index_river_wb
    res_index_f = network.river.reservoir_index_f
    lake_index_f = network.river.lake_index_f

    # total surface water demand for each land cell
    vertical.waterallocation.surfacewater_alloc .= 0.0
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
    # from reservoir and lake 
    if !isnothing(river.reservoir)
        @. river.abstraction[res_index_f] = 0.0
        @. river.reservoir.volume -=
            river.waterallocation.act_surfacewater_abst_vol[res_index_f]
    elseif !isnothing(river.lake)
        @. river.abstraction[lake_index_f] = 0.0
        @. river.lake.volume -=
            river.waterallocation.act_surfacewater_abst_vol[lake_index_f]
    end

    # local groundwater demand and allocation
    groundwater_allocation_local(vertical, lateral.subsurface, network.land.area)
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
    returnflow_livestock = livestock_return_flow(vertical)
    returnflow_domestic = domestic_return_flow(vertical)
    returnflow_industry = industry_return_flow(vertical)

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
