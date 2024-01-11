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
    nonirri_demand_gross::Vector{T}     # non irrigation gross demand
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

function initialize_water_allocation(n)
    waterallocation = WaterAllocation(
        irri_demand_gross = zeros(Float, n),
        nonirri_demand_gross = zeros(Float, n),
    )

    return waterallocation
end

function update_water_demand(sbm::SBM)
    for i = 1:sbm.n
        industry_dem = sbm.industry === nothing ? 0.0 : sbm.industry.demand_gross[i]
        domestic_dem = sbm.domestic === nothing ? 0.0 : sbm.domestic.demand_gross[i]
        livestock_dem = sbm.livestock === nothing ? 0.0 : sbm.livestock.demand_gross[i]

        nonirri_dem_gross = industry_dem + domestic_dem + livestock_dem

        irri_dem_gross = 0.0
        if sbm.nonpaddy !== nothing && sbm.nonpaddy.irrigation_areas[i] !== 0
            usl, _ = set_layerthickness(sbm.zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
            for k = 1:sbm.n_unsatlayers
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
        sbm.waterallocation.nonirri_dem_gross[i] = nonirri_dem_gross
    end
end
