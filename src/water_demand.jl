@get_units @with_kw struct Industry{T}
    demand_gross::Vector{T}                 # gross industry water demand [mm Δt⁻¹]
    demand_net::Vector{T}                   # net industry water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
end

@get_units @with_kw struct Domestic{T}
    demand_gross::Vector{T}                 # gross domestic water demand [mm Δt⁻¹]
    demand_net::Vector{T}                   # net domestic water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
end

@get_units @with_kw struct Livestock{T}
    demand_gross::Vector{T}                 # gross livestock water demand [mm Δt⁻¹]
    demand_net::Vector{T}                   # net livestock water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
end

@get_units @with_kw struct NonPaddy{T}
    demand_gross::Vector{T}                 # irrigation gross demand [mm Δt⁻¹] 
    irrigation_efficiency::Vector{T} | "-"  # irrigation efficiency [-]
    irrigation_areas::Vector{Bool} | "-"    # irrigation areas [-]
end

@get_units @with_kw struct Paddy{T}
    demand_gross::Vector{T}                 # irrigation gross demand [mm Δt⁻¹] 
    irrigation_efficiency::Vector{T} | "-"  # irrigation efficiency [-]
    irrigation_areas::Vector{Bool} | "-"    # irrigation areas [-]
    h_min::Vector{T} | "mm"                 # minimum required water depth in the irrigated rice field [mm]
    h_max::Vector{T} | "mm"                 # optimal water depth in the irrigated rice fields [mm]
    h_p::Vector{T} | "mm"                   # water depth when rice field starts spilling water (overflow) [mm]
    h::Vector{T} | "mm"                     # actual water depth in rice field [mm]
end

@get_units @with_kw struct WaterAllocation{T}
    irri_demand_gross::Vector{T}                        # irrigation gross demand [mm Δt⁻¹]
    nonirri_demand_gross::Vector{T}                     # non irrigation (industry and domestic) gross demand [mm Δt⁻¹]
    total_gross_demand::Vector{T}                       # total gross demand [mm Δt⁻¹]
    frac_sw_used::Vector{T} | "-"                       # fraction surface water used [-]
    areas::Vector{Int} | "-"                            # allocation areas [-]
    act_surfacewater_abst::Vector{T}                    # actual surface water abstraction [mm Δt⁻¹]
    act_surfacewater_abst_vol::Vector{T} | "m3 Δt-1"    # actual surface water abstraction [m³ Δt⁻¹]
    available_surfacewater::Vector{T} | "m3"            # available surface water [m³]
    surfacewater_demand::Vector{T}                      # demand from surface water [mm Δt⁻¹]
    surfacewater_alloc::Vector{T}                       # allocation from surface water [mm Δt⁻¹]
    act_groundwater_abst::Vector{T}                     # actual groundwater abstraction [mm Δt⁻¹]
    act_groundwater_abst_vol::Vector{T} | "m3 Δt-1"     # actual groundwater abstraction [m³ Δt⁻¹]
    available_groundwater::Vector{T} | "m3"             # available groundwater [m³]
    groundwater_demand::Vector{T}                       # demand from groundwater [mm Δt⁻¹]
    groundwater_alloc::Vector{T}                        # allocation from groundwater [mm Δt⁻¹]
    irri_alloc::Vector{T}                               # allocated water for irrigation [mm Δt⁻¹]
    nonirri_returnflow::Vector{T}                       # return flow from non irrigation [mm Δt⁻¹]
end

function set_returnflow_fraction(returnflow_fraction, demand_gross, demand_net)
    for i in eachindex(demand_gross)
        fraction = divide(demand_net[i], demand_gross[i])
        returnflow_fraction[i] = 1.0 - fraction
    end
    return returnflow_fraction
end

function initialize_domestic_demand(nc, config, inds, Δt)

    demand_gross =
        ncread(
            nc,
            config,
            "vertical.domestic.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.domestic.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)

    returnflow_fraction =
        set_returnflow_fraction(fill(mv, length(inds)), demand_gross, demand_net)

    domestic = Domestic{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
    )

    return domestic
end

function initialize_industry_demand(nc, config, inds, Δt)

    demand_gross =
        ncread(
            nc,
            config,
            "vertical.industry.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.industry.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)

    returnflow_fraction =
        set_returnflow_fraction(fill(mv, length(inds)), demand_gross, demand_net)

    industry = Industry{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
    )

    return industry
end

function initialize_livestock_demand(nc, config, inds, Δt)

    demand_gross =
        ncread(
            nc,
            config,
            "vertical.livestock.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.livestock.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)

    returnflow_fraction =
        set_returnflow_fraction(fill(mv, length(inds)), demand_gross, demand_net)

    livestock = Livestock{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
    )

    return livestock
end

function initialize_paddy(nc, config, inds)

    h_min =
        ncread(nc, config, "vertical.paddy.h_min"; sel = inds, defaults = 0.0, type = Float)
    h_max =
        ncread(nc, config, "vertical.paddy.h_max"; sel = inds, defaults = 0.0, type = Float)
    h_p = ncread(nc, config, "vertical.paddy.h_p"; sel = inds, defaults = 0.0, type = Float)
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
        defaults = 1,
        #optional = false,
        type = Int,
        #fill = 0,
    )

    paddy = Paddy{Float}(
        demand_gross = fill(mv, length(inds)),
        irrigation_efficiency = efficiency,
        h_min = h_min,
        h_max = h_max,
        h_p = h_p,
        irrigation_areas = areas,
        h = fill(0.0, length(inds)),
    )

    return paddy
end

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
        #optional = false,
        type = Int,
        #fill = 0,
    )

    nonpaddy = NonPaddy{Float}(
        demand_gross = fill(mv, length(inds)),
        irrigation_efficiency = efficiency,
        irrigation_areas = areas,
    )

    return nonpaddy
end

function initialize_water_allocation(nc, config, inds, nriv)

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

    waterallocation = WaterAllocation(
        irri_demand_gross = zeros(Float, n),
        nonirri_demand_gross = zeros(Float, n),
        total_gross_demand = zeros(Float, n),
        frac_sw_used = frac_sw_used,
        areas = areas,
        act_surfacewater_abst = zeros(Float, nriv),
        act_surfacewater_abst_vol = zeros(Float, nriv),
        available_surfacewater = zeros(Float, nriv),
        surfacewater_demand = zeros(Float, n),
        surfacewater_alloc = zeros(Float, n),
        act_groundwater_abst = zeros(Float, n),
        act_groundwater_abst_vol = zeros(Float, n),
        available_groundwater = zeros(Float, n),
        groundwater_demand = zeros(Float, n),
        groundwater_alloc = zeros(Float, n),
        irri_alloc = zeros(Float, n),
        nonirri_returnflow = zeros(Float, n),
    )

    return waterallocation
end

function update_water_demand(sbm::SBM)
    for i = 1:sbm.n
        industry_dem = sbm.industry === nothing ? 0.0 : sbm.industry.demand_gross[i]
        domestic_dem = sbm.domestic === nothing ? 0.0 : sbm.domestic.demand_gross[i]
        livestock_dem = sbm.livestock === nothing ? 0.0 : sbm.livestock.demand_gross[i]

        irri_dem_gross = 0.0
        if sbm.nonpaddy !== nothing && sbm.nonpaddy.irrigation_areas[i] !== 0
            usl, _ = set_layerthickness(sbm.zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
            # TODO: only include root zone
            for k = 1:sbm.n_unsatlayers[i]
                vwc_fc =
                    vwc_brooks_corey(-100.0, sbm.hb[i], sbm.θₛ[i], sbm.θᵣ[i], sbm.c[i][k])
                vwc_h3 = vwc_brooks_corey(
                    sbm.h3[i],
                    sbm.hb[i],
                    sbm.θₛ[i],
                    sbm.θᵣ[i],
                    sbm.c[i][k],
                )
                depletion = (vwc_fc * usl[k]) - sbm.ustorelayerdepth[i][k]
                raw = (vwc_fc - vwc_h3) * usl[k] # readily available water
                if depletion >= raw
                    irri_dem_gross += depletion
                end
                # limit irrigation demand to infiltration capacity    
                infiltration_capacity =
                    sbm.soilinfredu[i] * (sbm.infiltcappath[i] + sbm.infiltcapsoil[i])
                irri_dem_gross = min(irri_dem_gross, infiltration_capacity)
                irri_dem_gross /= sbm.nonpaddy.irrigation_efficiency[i]
            end
        elseif sbm.paddy !== nothing && sbm.paddy.irrigation_areas[i] !== 0
            irr_depth_paddy =
                sbm.paddy.h[i] < sbm.paddy.h_min[i] ?
                (sbm.paddy.h_max[i] - sbm.paddy.h[i]) : 0.0
            irri_dem_gross += irr_depth_paddy / sbm.paddy.irrigation_efficiency[i]
        end
        sbm.waterallocation.irri_demand_gross[i] = irri_dem_gross
        sbm.waterallocation.nonirri_demand_gross[i] =
            industry_dem + domestic_dem + livestock_dem
        sbm.waterallocation.total_gross_demand[i] =
            irri_dem_gross + industry_dem + domestic_dem + livestock_dem
    end
end

function update_water_allocation(model)

    @unpack network, lateral, vertical = model
    @unpack waterallocation = vertical

    inds_river = network.river.indices_allocation_areas
    inds_land = network.land.indices_allocation_areas

    index_river = network.land.index_river

    n = length(network.land.indices)
    # local surface water abstraction (river, excluding reservoirs and lakes)
    for i = 1:n
        waterallocation.surfacewater_alloc[i] = 0.0
        waterallocation.surfacewater_demand[i] =
            waterallocation.frac_sw_used[i] * waterallocation.nonirri_demand_gross[i] +
            waterallocation.frac_sw_used[i] * waterallocation.irri_demand_gross[i]

        if index_river[i] > 0.0 # TODO: exclude reservoir/lake cells
            # check for abstraction through inflow and adjust available volume
            if lateral.river.inflow[index_river[i]] < 0.0
                inflow = lateral.river.inflow[index_river[i]] * vertical.Δt
                available_volume =
                    max(lateral.river.volume[index_river[i]] * 0.80 + inflow, 0.0)
            else
                available_volume = lateral.river.volume[index_river[i]] * 0.80
            end
            surfacewater_demand_vol =
                waterallocation.surfacewater_demand[i] * 0.001 * network.land.area[i]
            abstraction_vol = min(surfacewater_demand_vol, available_volume)
            waterallocation.available_surfacewater[index_river[i]] =
                available_volume - abstraction_vol
            abstraction = (abstraction_vol / network.land.area[i]) * 1000.0
            waterallocation.surfacewater_demand[i] -= abstraction
            waterallocation.act_surfacewater_abst[index_river[i]] = abstraction
            waterallocation.surfacewater_alloc[i] = abstraction
        end
    end

    # surface water demand and availability for allocation areas
    m = length(inds_river)
    for i = 1:m
        # surface water demand (allocation area)
        sw_demand_vol = 0.0
        for j in inds_land[i]
            sw_demand_vol +=
                waterallocation.surfacewater_demand[j] * 0.001 * network.land.area[j]
        end
        # surface water availability (allocation area)
        sw_available = 0.0
        for j in inds_river[i]
            if lateral.river.reservoir_index[j] > 0
                k = lateral.river.reservoir_index[j]
                waterallocation.available_surfacewater[j] =
                    lateral.river.reservoir.volume[k] * 0.98
                sw_available += waterallocation.available_surfacewater[j]
            elseif lateral.river.lake_index[j] > 0
                k = lateral.river.lake_index[j]
                waterallocation.available_surfacewater[j] =
                    lateral.river.lake.storage[k] * 0.98
                sw_available += waterallocation.available_surfacewater[j]
            else
                sw_available += waterallocation.available_surfacewater[j]
            end
        end
        sw_abstraction = min(sw_available, sw_demand_vol)

        frac_abstract_sw = divide(sw_abstraction, sw_available)
        frac_allocate_sw = divide(sw_abstraction, sw_demand_vol)

        for j in inds_river[i]
            waterallocation.act_surfacewater_abst_vol[j] +=
                frac_abstract_sw * waterallocation.available_surfacewater[j]
            waterallocation.act_surfacewater_abst[j] =
                waterallocation.act_surfacewater_abst_vol[j] / network.river.area[j]
        end

        for j in inds_land[i]
            waterallocation.surfacewater_alloc[j] +=
                frac_allocate_sw * waterallocation.surfacewater_demand[j]
        end
    end
    @. lateral.river.abstraction = waterallocation.act_surfacewater_abst_vol / vertical.Δt

    # local groundwater abstraction
    for i = 1:n
        waterallocation.groundwater_demand[i] =
            waterallocation.irri_demand_gross[i] + waterallocation.nonirri_demand_gross[i] -
            waterallocation.surfacewater_alloc[i]
        groundwater_demand_vol =
            waterallocation.groundwater_demand[i] * 0.001 * network.land.area[i]
        available_volume = lateral.subsurface.volume[i] * 0.75
        abstraction_vol = min(groundwater_demand_vol, available_volume)
        waterallocation.available_groundwater[i] = available_volume - abstraction_vol
        abstraction = (abstraction_vol / network.land.area[i]) * 1000.0
        waterallocation.groundwater_demand[i] -= abstraction
        waterallocation.act_groundwater_abst[i] = abstraction
        waterallocation.groundwater_alloc[i] = abstraction
    end
    # groundwater demand and availability for allocation areas
    for i = 1:m
        gw_demand_vol = 0.0
        gw_available = 0.0
        for j in inds_land[i]
            gw_demand_vol +=
                waterallocation.groundwater_demand[j] * 0.001 * network.land.area[j]
            gw_available += waterallocation.available_groundwater[j]
        end
        gw_abstraction = min(gw_available, gw_demand_vol)

        frac_abstract_gw = divide(gw_abstraction, gw_available)
        frac_allocate_gw = divide(gw_abstraction, gw_demand_vol)

        for j in inds_land[i]
            waterallocation.act_groundwater_abst_vol[j] +=
                frac_abstract_gw * waterallocation.available_groundwater[j]
            waterallocation.act_groundwater_abst[j] =
                waterallocation.act_groundwater_abst[j] / network.land.area[j]
            waterallocation.groundwater_alloc[j] +=
                frac_allocate_gw * waterallocation.groundwater_demand[j]
        end
    end

    for i = 1:n
        total_alloc =
            waterallocation.groundwater_alloc[i] + waterallocation.surfacewater_alloc[i]
        frac_irri = divide(
            waterallocation.irri_demand_gross[i],
            waterallocation.total_gross_demand[i],
        )
        waterallocation.irri_alloc[i] = frac_irri * total_alloc

        nonirri_alloc = total_alloc - waterallocation.irri_alloc[i]
        frac_livestock = divide(
            vertical.livestock.demand_gross[i],
            waterallocation.nonirri_demand_gross[i],
        )
        livestock_alloc = frac_livestock * nonirri_alloc
        frac_domestic = divide(
            vertical.domestic.demand_gross[i],
            waterallocation.nonirri_demand_gross[i],
        )
        domestic_alloc = frac_domestic * nonirri_alloc
        frac_industry = divide(
            vertical.industry.demand_gross[i],
            waterallocation.nonirri_demand_gross[i],
        )
        industry_alloc = frac_industry * nonirri_alloc

        waterallocation.nonirri_returnflow[i] =
            vertical.livestock.returnflow_fraction[i] * livestock_alloc +
            vertical.domestic.returnflow_fraction[i] * domestic_alloc +
            vertical.industry.returnflow_fraction[i] * industry_alloc
    end
end
