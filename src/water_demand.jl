@get_units @with_kw struct Industry{T}
    demand_gross::Vector{T}   # gross industry water demand [mm Δt⁻¹]
    demand_net::Vector{T}     # net industry water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

@get_units @with_kw struct Domestic{T}
    demand_gross::Vector{T}   # gross domestic water demand [mm Δt⁻¹]
    demand_net::Vector{T}     # net domestic water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

@get_units @with_kw struct Livestock{T}
    demand_gross::Vector{T}  # gross livestock water demand [mm Δt⁻¹]
    demand_net::Vector{T}    # net livestock water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

@get_units @with_kw struct NonPaddy{T}
    demand_gross::Vector{T}
    irrigation_efficiency::Vector{T}
    irrigation_areas::Vector{Bool}
end

@get_units @with_kw struct Paddy{T}
    demand_gross::Vector{T}
    irrigation_efficiency::Vector{T}
    irrigation_areas::Vector{Bool}
    h_min::Vector{T}                    # minimum required water depth in the irrigated rice field [mm]
    h_max::Vector{T}                    # optimal water depth in the irrigated rice fields [mm]
    h_p::Vector{T}                      # water depth when rice field starts spilling water (overflow) [mm]
    h::Vector{T}                        # actual water depth in rice field
end

@get_units @with_kw struct WaterAllocation{T}
    irri_demand_gross::Vector{T}        # irrigation gross demand
    nonirri_demand_gross::Vector{T}     # non irrigation (industry and domestic) gross demand
    total_gross_demand::Vector{T}       # total gross demand
    frac_sw_used::Vector{T}             # fraction surface water used
    areas::Vector{Int}                  # allocation areas
    act_surfacewater_abst::Vector{T}    # actual surface water abstraction
    available_surfacewater::Vector{T}   # available surface water
    surfacewater_demand::Vector{T}      # demand from surface water
    surfacewater_alloc::Vector{T}       # allocation from surface water
    act_groundwater_abst::Vector{T}     # actual groundwater abstraction
    available_groundwater::Vector{T}    # available groundwater
    groundwater_demand::Vector{T}       # demand from groundwater
    groundwater_alloc::Vector{T}        # allocation from groundwater
end

function set_returnflow_fraction(returnflow_fraction, demand_gross, demand_net)
    for i in eachindex(demand_gross)
        fraction = demand_net[i] / demand_gross[i]
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
        available_surfacewater = zeros(Float, nriv),
        surfacewater_demand = zeros(Float, n),
        surfacewater_alloc = zeros(Float, n),
        act_groundwater_abst = zeros(Float, n),
        available_groundwater = zeros(Float, n),
        groundwater_demand = zeros(Float, n),
        groundwater_alloc = zeros(Float, n),
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
            0.001 *
            (
                waterallocation.frac_sw_used[i] * waterallocation.nonirri_demand_gross[i] +
                waterallocation.frac_sw_used[i] * waterallocation.irri_demand_gross[i]
            ) *
            network.land.area[i]
        if index_river[i] > 0.0 # TODO: exclude reservoir/lake cells
            # check for abstraction through inflow and adjust available volume
            if lateral.river.inflow[index_river[i]] < 0.0
                inflow = lateral.river.inflow[index_river[i]] * vertical.Δt
                available_volume =
                    max(lateral.river.volume[index_river[i]] * 0.80 + inflow, 0.0)
            else
                available_volume = lateral.river.volume[index_river[i]] * 0.80
            end
            abstraction = min(waterallocation.surfacewater_demand[i], available_volume)
            waterallocation.available_surfacewater[index_river[i]] =
                available_volume - abstraction
            waterallocation.surfacewater_demand[i] -= abstraction
            waterallocation.act_surfacewater_abst[index_river[i]] = abstraction
            waterallocation.surfacewater_alloc[i] = abstraction
        end
    end

    # surface water demand and availability for allocation areas
    m = length(inds_river)
    for i = 1:m
        # surface water demand (allocation area)
        sw_demand = 0.0
        for j in inds_land[i]
            sw_demand += waterallocation.surfacewater_demand[j]
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
        sw_abstraction = min(sw_available, sw_demand)

        frac_abstract_sw =
            sw_available > 0.0 ? min(sw_abstraction / sw_available, 1.0) : 0.0
        frac_allocate_sw = sw_demand > 0.0 ? min(sw_abstraction / sw_demand, 1.0) : 0.0

        for j in inds_river[i]
            waterallocation.act_surfacewater_abst[j] +=
                frac_abstract_sw * waterallocation.available_surfacewater[j]
        end

        for j in inds_land[i]
            waterallocation.surfacewater_alloc[j] +=
                frac_allocate_sw * waterallocation.surfacewater_demand[j]
        end
    end
    @. lateral.river.abstraction = waterallocation.act_surfacewater_abst / vertical.Δt

    # local groundwater abstraction
    for i = 1:n
        waterallocation.groundwater_demand[i] =
            0.001 *
            (
                waterallocation.irri_demand_gross[i] +
                waterallocation.nonirri_demand_gross[i]
            ) *
            network.land.area[i] - waterallocation.surfacewater_alloc[i]

        available_volume = lateral.subsurface.volume[i] * 0.75
        abstraction = min(waterallocation.groundwater_demand[i], available_volume)
        waterallocation.available_groundwater[i] = available_volume - abstraction
        waterallocation.groundwater_demand[i] -= abstraction
        waterallocation.act_groundwater_abst[i] = abstraction
        waterallocation.groundwater_alloc[i] = abstraction
    end
    # groundwater demand and availability for allocation areas
    for i = 1:m
        gw_demand = 0.0
        gw_available = 0.0
        for j in inds_land[i]
            gw_demand += waterallocation.groundwater_demand[j]
            gw_available += waterallocation.available_groundwater[j]
        end
        gw_abstraction = min(gw_available, gw_demand)

        frac_abstract_gw =
            gw_available > 0.0 ? min(gw_abstraction / gw_available, 1.0) : 0.0
        frac_allocate_gw = gw_demand > 0.0 ? min(gw_abstraction / gw_demand, 1.0) : 0.0

        for j in inds_land[i]
            waterallocation.act_groundwater_abst[j] +=
                frac_abstract_gw * waterallocation.available_groundwater[j]
            waterallocation.groundwater_alloc[j] +=
                frac_allocate_gw * waterallocation.groundwater_demand[j]
        end
    end
end
