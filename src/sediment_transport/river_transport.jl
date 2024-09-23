abstract type AbstractSedimentRiverTransportModel end

## Total sediment transport in overland flow structs and functions
@get_units @with_kw struct SedimentRiverTransportVars{T}
    # Sediment flux [ton]
    amount::Vector{T} | "t dt-1"
    clay::Vector{T} | "t dt-1"
    silt::Vector{T} | "t dt-1"
    sand::Vector{T} | "t dt-1"
    sagg::Vector{T} | "t dt-1"
    lagg::Vector{T} | "t dt-1"
    gravel::Vector{T} | "t dt-1"
    # Total Sediment deposition [ton]
    deposition::Vector{T} | "t dt-1"
    # Total sediment erosion (from store + direct river bed/bank) [ton]
    erosion::Vector{T} | "t dt-1"
    # Sediment / particle left in the cell [ton] - states
    leftover_clay::Vector{T} | "t dt-1"
    leftover_silt::Vector{T} | "t dt-1"
    leftover_sand::Vector{T} | "t dt-1"
    leftover_sagg::Vector{T} | "t dt-1"
    leftover_lagg::Vector{T} | "t dt-1"
    leftover_gravel::Vector{T} | "t dt-1"
    # Sediment / particle stored on the river bed after deposition [ton] -states
    store_clay::Vector{T} | "t dt-1"
    store_silt::Vector{T} | "t dt-1"
    store_sand::Vector{T} | "t dt-1"
    store_sagg::Vector{T} | "t dt-1"
    store_lagg::Vector{T} | "t dt-1"
    store_gravel::Vector{T} | "t dt-1"
end

function sediment_river_transport_vars(n)
    vars = SedimentRiverTransportVars(;
        amount = fill(mv, n),
        clay = fill(mv, n),
        silt = fill(mv, n),
        sand = fill(mv, n),
        sagg = fill(mv, n),
        lagg = fill(mv, n),
        gravel = fill(mv, n),
        deposition = fill(mv, n),
        erosion = fill(mv, n),
        leftover_clay = fill(mv, n),
        leftover_silt = fill(mv, n),
        leftover_sand = fill(mv, n),
        leftover_sagg = fill(mv, n),
        leftover_lagg = fill(mv, n),
        leftover_gravel = fill(mv, n),
        store_clay = fill(mv, n),
        store_silt = fill(mv, n),
        store_sand = fill(mv, n),
        store_sagg = fill(mv, n),
        store_lagg = fill(mv, n),
        store_gravel = fill(mv, n),
    )
    return vars
end

@get_units @with_kw struct SedimentRiverTransportBC{T}
    # Waterlevel
    waterlevel::Vector{T} | "t dt-1"
    # Discharge
    q::Vector{T} | "m3 s-1"
    # Transport capacity of the flow
    transport_capacity::Vector{T} | "t dt-1"
    # Sediment input from land erosion
    erosion_land_clay::Vector{T} | "t dt-1"
    erosion_land_silt::Vector{T} | "t dt-1"
    erosion_land_sand::Vector{T} | "t dt-1"
    erosion_land_sagg::Vector{T} | "t dt-1"
    erosion_land_lagg::Vector{T} | "t dt-1"
    # Sediment available from direct river erosion
    potential_erosion_river::Vector{T} | "t dt-1"
end

function sediment_river_transport_bc(n)
    bc = SedimentRiverTransportBC(;
        waterlevel = fill(mv, n),
        q = fill(mv, n),
        transport_capacity = fill(mv, n),
        erosion_land_clay = fill(mv, n),
        erosion_land_silt = fill(mv, n),
        erosion_land_sand = fill(mv, n),
        erosion_land_sagg = fill(mv, n),
        erosion_land_lagg = fill(mv, n),
        potential_erosion_river = fill(mv, n),
    )
    return bc
end

# Parameters for river transport
@get_units @with_kw struct SedimentRiverTransportParameters{T}
    # River bed/bank content clay
    clay_fraction::Vector{T} | "-"
    # River bed/bank content silt
    silt_fraction::Vector{T} | "-"
    # River bed/bank content sand
    sand_fraction::Vector{T} | "-"
    # River bed/bank content gravel
    gravel_fraction::Vector{T} | "-"
    # Clay mean diameter
    dm_clay::Vector{T} | "µm"
    # Silt mean diameter
    dm_silt::Vector{T} | "µm"
    # Sand mean diameter
    dm_sand::Vector{T} | "µm"
    # Small aggregates mean diameter
    dm_sagg::Vector{T} | "µm"
    # Large aggregates mean diameter
    dm_lagg::Vector{T} | "µm"
    # Gravel mean diameter
    dm_gravel::Vector{T} | "µm"
    # Waterbodies outlets
    waterbodies_locs::Vector{Bool} | "-"
    # Waterbodies area
    waterbodies_area::Vector{T} | "m2"
    # Waterbodies trapping efficiency
    waterbodies_trapping_efficiency::Vector{T} | "-"
end

function initialize_sediment_river_transport_params(nc, config, inds)
    n = length(inds)
    clay_fraction = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.clay_fraction";
        sel = inds,
        defaults = 0.15,
        type = Float,
    )
    silt_fraction = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.silt_fraction";
        sel = inds,
        defaults = 0.65,
        type = Float,
    )
    sand_fraction = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.sand_fraction";
        sel = inds,
        defaults = 0.15,
        type = Float,
    )
    gravel_fraction = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.sagg_fraction";
        sel = inds,
        defaults = 0.05,
        type = Float,
    )
    # Check that river fractions sum to 1
    river_fractions = clay_fraction + silt_fraction + sand_fraction + gravel_fraction
    if any(abs.(river_fractions .- 1.0) .> 1e-3)
        error("Particle fractions in the river bed must sum to 1")
    end
    dm_clay = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.dm_clay";
        sel = inds,
        defaults = 2.0,
        type = Float,
    )
    dm_silt = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.dm_silt";
        sel = inds,
        defaults = 10.0,
        type = Float,
    )
    dm_sand = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.dm_sand";
        sel = inds,
        defaults = 200.0,
        type = Float,
    )
    dm_sagg = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.dm_sagg";
        sel = inds,
        defaults = 30.0,
        type = Float,
    )
    dm_lagg = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.dm_lagg";
        sel = inds,
        defaults = 500.0,
        type = Float,
    )
    dm_gravel = ncread(
        nc,
        config,
        "lateral.river.sediment_flux.parameters.dm_gravel";
        sel = inds,
        defaults = 2000.0,
        type = Float,
    )
    # Waterbodies
    wblocs = zeros(Float, n)
    wbarea = zeros(Float, n)
    wbtrap = zeros(Float, n)
    do_reservoirs = get(config.model, "doreservoir", false)::Bool
    do_lakes = get(config.model, "dolake", false)::Bool

    if do_reservoirs
        reslocs = ncread(
            nc,
            config,
            "lateral.river.sediment_flux.parameters.reslocs";
            optional = false,
            sel = inds,
            type = Float,
            fill = 0,
        )
        resarea = ncread(
            nc,
            config,
            "lateral.river.sediment_flux.parameters.resarea";
            optional = false,
            sel = inds,
            type = Float,
            fill = 0.0,
        )
        restrapefficiency = ncread(
            nc,
            config,
            "lateral.river.sediment_flux.parameters.restrapeff";
            optional = false,
            sel = inds,
            type = Float,
            defaults = 1.0,
            fill = 0.0,
        )
        wblocs = wblocs .+ reslocs
        wbarea = wbarea .+ resarea
        wbtrap = wbtrap .+ restrapefficiency
    end

    if do_lakes
        lakelocs = ncread(
            nc,
            config,
            "lateral.river.sediment_flux.parameters.lakelocs";
            optional = false,
            sel = inds,
            type = Float,
            fill = 0,
        )
        lakearea = ncread(
            nc,
            config,
            "lateral.river.sediment_flux.parameters.lakearea";
            optional = false,
            sel = inds,
            type = Float,
            fill = 0.0,
        )
        wblocs = wblocs .+ lakelocs
        wbarea = wbarea .+ lakearea
    end

    river_parameters = SedimentRiverTransportParameters(;
        clay_fraction = clay_fraction,
        silt_fraction = silt_fraction,
        sand_fraction = sand_fraction,
        gravel_fraction = gravel_fraction,
        dm_clay = dm_clay,
        dm_silt = dm_silt,
        dm_sand = dm_sand,
        dm_sagg = dm_sagg,
        dm_lagg = dm_lagg,
        dm_gravel = dm_gravel,
        waterbodies_locs = wblocs .> 0,
        waterbodies_area = wbarea,
        waterbodies_trapping_efficiency = wbtrap,
    )

    return river_parameters
end

@get_units @with_kw struct SedimentRiverTransportModel{T} <:
                           AbstractSedimentRiverTransportModel
    boundary_conditions::SedimentRiverTransportBC{T} | "-"
    parameters::SedimentRiverTransportParameters{T} | "-"
    variables::SedimentRiverTransportVars{T} | "-"
end

function initialize_sediment_river_transport_model(nc, config, inds)
    n = length(inds)
    vars = sediment_river_transport_vars(n)
    params = initialize_sediment_river_transport_params(nc, config, inds)
    bc = sediment_river_transport_bc(n)
    model = SedimentRiverTransportModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::SedimentRiverTransportModel, network, geometry::RiverGeometry, ts)
    (;
        waterlevel,
        q,
        transport_capacity,
        erosion_land_clay,
        erosion_land_silt,
        erosion_land_sand,
        erosion_land_sagg,
        erosion_land_lagg,
        potential_erosion_river_bed,
        potential_erosion_river_bank,
    ) = model.boundary_conditions
    (; clay_fraction, silt_fraction, sand_fraction, gravel_fraction) = model.parameters
    (;
        amount,
        clay,
        silt,
        sand,
        sagg,
        lagg,
        gravel,
        deposition,
        erosion,
        leftover_clay,
        leftover_silt,
        leftover_sand,
        leftover_sagg,
        leftover_lagg,
        leftover_gravel,
        store_clay,
        store_silt,
        store_sand,
        store_sagg,
        store_lagg,
        store_gravel,
    ) = model.variables

    @unpack graph, order = network

    # Sediment transport - water balance in the river
    for v in order
        ### Sediment input in the cell (left from previous timestep + from land + from upstream outflux) ###
        input_clay = leftover_clay[v] + erosion_land_clay[v]
        input_silt = leftover_silt[v] + erosion_land_silt[v]
        input_sand = leftover_sand[v] + erosion_land_sand[v]
        input_sagg = leftover_sagg[v] + erosion_land_sagg[v]
        input_lagg = leftover_lagg[v] + erosion_land_lagg[v]
        input_gravel = leftover_gravel[v]

        # Add upstream contribution
        upstream_nodes = inneighbors(graph, v)
        if !isempty(upstream_nodes)
            for i in upstream_nodes
                if clay[i] >= 0.0 # avoid NaN from upstream non-river cells
                    input_clay += clay[i]
                    input_silt += silt[i]
                    input_sand += sand[i]
                    input_sagg += sagg[i]
                    input_lagg += lagg[i]
                    input_gravel += gravel[i]
                end
            end
        end

        input_sediment =
            input_clay + input_silt + input_sand + input_sagg + input_lagg + input_gravel

        ### River erosion ###
        # Erosion only if the load is below the transport capacity of the flow.
        sediment_need = max(transport_capacity - input_sediment, 0.0)
        # No erosion in reservoirs
        if waterbodies[v]
            sediment_need = 0.0
        end

        # Available sediment stored from previous deposition
        store_sediment =
            store_clay[v] +
            store_silt[v] +
            store_sand[v] +
            store_sagg[v] +
            store_lagg[v] +
            store_gravel[v]

        # Direct erosion from the river bed/bank
        if sediment_need > store_sediment
            # Effective sediment needed fom river bed and bank erosion [ton]
            effsediment_need = sediment_need - store_sediment
            # Relative potential erosion rates of the bed and the bank [-]
            if (potential_erosion_river_bank + potential_erosion_river_bed > 0.0)
                RTEbank =
                    potential_erosion_river_bank /
                    (potential_erosion_river_bank + potential_erosion_river_bed)
            else
                RTEbank = 0.0
            end
            RTEbed = 1.0 - RTEbank

            # Actual bed and bank erosion
            erosion_bank = max(RTEbank * effsediment_need, potential_erosion_river_bank)
            erosion_bed = max(RTEbed * effsediment_need, potential_erosion_river_bed)
            erosion_river = erosion_bank + erosion_bed
            # Per particle
            erosion_clay = erosion_river * clay_fraction[v]
            erosion_silt = erosion_river * silt_fraction[v]
            erosion_sand = erosion_river * sand_fraction[v]
            erosion_gravel = erosion_river * gravel_fraction[v]
            # No small and large aggregates in the river bed/bank
            erosion_sagg = 0.0
            erosion_lagg = 0.0
        else
            erosion_clay = 0.0
            erosion_silt = 0.0
            erosion_sand = 0.0
            erosion_gravel = 0.0
            erosion_sagg = 0.0
            erosion_lagg = 0.0
        end

        # Erosion/degradation of the previously deposited sediment (from clay to gravel) [ton]
        if sediment_need > 0.0
            # Erosion in priority of the smaller particles
            # Clay
            if store_clay[v] > 0.0
                erosion_store_clay, sediment_need, store_clay[v] =
                    river_erosion_store(sediment_need, store_clay[v])
                # Update the clay erosion
                erosion_clay += erosion_store_clay
            end
            # Silt
            if store_silt[v] > 0.0
                erosion_store_silt, sediment_need, store_silt[v] =
                    river_erosion_store(sediment_need, store_silt[v])
                # Update the silt erosion
                erosion_silt += erosion_store_silt
            end
            # Small aggregates
            if store_sagg[v] > 0.0
                erosion_store_sagg, sediment_need, store_sagg[v] =
                    river_erosion_store(sediment_need, store_sagg[v])
                # Update the sagg erosion
                erosion_sagg += erosion_store_sagg
            end
            # Sand
            if store_sand[v] > 0.0
                erosion_store_sand, sediment_need, store_sand[v] =
                    river_erosion_store(sediment_need, store_sand[v])
                # Update the sand erosion
                erosion_sand += erosion_store_sand
            end
            # Large aggregates
            if store_lagg[v] > 0.0
                erosion_store_lagg, sediment_need, store_lagg[v] =
                    river_erosion_store(sediment_need, store_lagg[v])
                # Update the lagg erosion
                erosion_lagg += erosion_store_lagg
            end
            # Gravel
            if store_gravel[v] > 0.0
                erosion_store_gravel, sediment_need, store_gravel[v] =
                    river_erosion_store(sediment_need, store_gravel[v])
                # Update the gravel erosion
                erosion_gravel += erosion_store_gravel
            end
        end

        # Compute total erosion
        erosion[v] =
            erosion_clay +
            erosion_silt +
            erosion_sand +
            erosion_sagg +
            erosion_lagg +
            erosion_gravel

        ### Deposition / settling ###
        # Different deposition if waterbody outlet or river
        if waterbodies_locs[v]
            # Deposition in waterbodies outlets
            deposition_clay = reservoir_deposition_camp(
                (input_clay + erosion_clay),
                q[v],
                waterlevel[v],
                waterbodies_area[v],
                waterbodies_trapping_efficiency[v],
                dm_clay[v],
                geometry.slope[v],
            )
            deposition_silt = reservoir_deposition_camp(
                (input_silt + erosion_silt),
                q[v],
                waterlevel[v],
                waterbodies_area[v],
                waterbodies_trapping_efficiency[v],
                dm_silt[v],
                geometry.slope[v],
            )
            deposition_sand = reservoir_deposition_camp(
                (input_sand + erosion_sand),
                q[v],
                waterlevel[v],
                waterbodies_area[v],
                waterbodies_trapping_efficiency[v],
                dm_sand[v],
                geometry.slope[v],
            )
            deposition_sagg = reservoir_deposition_camp(
                (input_sagg + erosion_sagg),
                q[v],
                waterlevel[v],
                waterbodies_area[v],
                waterbodies_trapping_efficiency[v],
                dm_sagg[v],
                geometry.slope[v],
            )
            deposition_lagg = reservoir_deposition_camp(
                (input_lagg + erosion_lagg),
                q[v],
                waterlevel[v],
                waterbodies_area[v],
                waterbodies_trapping_efficiency[v],
                dm_lagg[v],
                geometry.slope[v],
            )
            deposition_gravel = reservoir_deposition_camp(
                (input_gravel + erosion_gravel),
                q[v],
                waterlevel[v],
                waterbodies_area[v],
                waterbodies_trapping_efficiency[v],
                dm_gravel[v],
                geometry.slope[v],
            )
        elseif waterbodies[v]
            # No deposition in waterbodies, only at the outlets
            deposition_clay = 0.0
            deposition_silt = 0.0
            deposition_sand = 0.0
            deposition_sagg = 0.0
            deposition_lagg = 0.0
            deposition_gravel = 0.0
        else
            # Deposition in the river
            # From trasnport capacity exceedance
            excess_sediment = max(input_sediment - transport_capacity, 0.0)
            if excess_sediment > 0.0
                # Sediment deposited in the channel (from gravel to clay) [ton]
                # Gravel
                deposition_gravel = ifelse(
                    excess_sediment > (input_gravel + erosion_gravel),
                    (input_gravel + erosion_gravel),
                    excess_sediment,
                )
                excess_sediment = max(excess_sediment - deposition_gravel, 0.0)
                # Large aggregates
                deposition_lagg = ifelse(
                    excess_sediment > (input_lagg + erosion_lagg),
                    (input_lagg + erosion_lagg),
                    excess_sediment,
                )
                excess_sediment = max(excess_sediment - deposition_lagg, 0.0)
                # Sand
                deposition_sand = ifelse(
                    excess_sediment > (input_sand + erosion_sand),
                    (input_sand + erosion_sand),
                    excess_sediment,
                )
                excess_sediment = max(excess_sediment - deposition_sand, 0.0)
                # Small aggregates
                deposition_sagg = ifelse(
                    excess_sediment > (input_sagg + erosion_sagg),
                    (input_sagg + erosion_sagg),
                    excess_sediment,
                )
                excess_sediment = max(excess_sediment - deposition_sagg, 0.0)
                # Silt
                deposition_silt = ifelse(
                    excess_sediment > (input_silt + erosion_silt),
                    (input_silt + erosion_silt),
                    excess_sediment,
                )
                excess_sediment = max(excess_sediment - deposition_silt, 0.0)
                # Clay
                deposition_clay = ifelse(
                    excess_sediment > (input_clay + erosion_clay),
                    (input_clay + erosion_clay),
                    excess_sediment,
                )
                excess_sediment = max(excess_sediment - deposition_clay, 0.0)
            else
                # Natural deposition from Einstein's formula (density controlled)
                # Particle fall velocity [m/s] from Stokes
                xs = ifelse(
                    q[v] > 0.0,
                    1.055 * geometry.length[v] / (q[v] / geometry.width[v]),
                    0.0,
                )
                xclay =
                    min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (dm_clay[v] / 1000)^2 / 3600)))
                deposition_clay = xclay * (input_clay + erosion_clay)
                xsilt =
                    min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (dm_silt[v] / 1000)^2 / 3600)))
                deposition_silt = xsilt * (input_silt + erosion_silt)
                xsand =
                    min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (dm_sand[v] / 1000)^2 / 3600)))
                deposition_sand = xsand * (input_sand + erosion_sand)
                xsagg =
                    min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (dm_sagg[v] / 1000)^2 / 3600)))
                deposition_sagg = xsagg * (input_sagg + erosion_sagg)
                xlagg =
                    min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (dm_lagg[v] / 1000)^2 / 3600)))
                deposition_lagg = xlagg * (input_lagg + erosion_lagg)
                xgrav =
                    min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (dm_grav[v] / 1000)^2 / 3600)))
                deposition_gravel = xgrav * (input_gravel + erosion_gravel)
            end
        end

        # Update the sediment store
        store_clay[v] += deposition_clay
        store_silt[v] += deposition_silt
        store_sand[v] += deposition_sand
        store_sagg[v] += deposition_sagg
        store_lagg[v] += deposition_lagg
        store_gravel[v] += deposition_gravel

        # Compute total deposition
        deposition[v] =
            deposition_clay +
            deposition_silt +
            deposition_sand +
            deposition_sagg +
            deposition_lagg +
            deposition_gravel

        ### Output loads ###
        # Sediment transported out of the cell during the timestep [ton]
        # 0 in case all sediment are deposited in the cell
        # Reduce the fraction so that there is still some sediment staying in the river cell
        fwaterout =
            min(q[v] * ts / (waterlevel[v] * geometry.width[v] * geometry.length[v]), 1.0)
        clay = fwaterout * (input_clay + erosion_clay - deposition_clay)
        silt = fwaterout * (input_silt + erosion_silt - deposition_silt)
        sand = fwaterout * (input_sand + erosion_sand - deposition_sand)
        sagg = fwaterout * (input_sagg + erosion_sagg - deposition_sagg)
        lagg = fwaterout * (input_lagg + erosion_lagg - deposition_lagg)
        gravel = fwaterout * (input_gravel + erosion_gravel - deposition_gravel)

        amount = clay + silt + sand + sagg + lagg + gravel

        ### Leftover / mass balance ###
        # Sediment left in the cell [ton]
        leftover_clay[v] = input_clay + erosion_clay - deposition_clay - clay
        leftover_silt[v] = input_silt + erosion_silt - deposition_silt - silt
        leftover_sand[v] = input_sand + erosion_sand - deposition_sand - sand
        leftover_sagg[v] = input_sagg + erosion_sagg - deposition_sagg - sagg
        leftover_lagg[v] = input_lagg + erosion_lagg - deposition_lagg - lagg
        leftover_gravel[v] = input_gravel + erosion_gravel - deposition_gravel - gravel
    end
end

abstract type AbstractSedimentConcentrationsRiverModel end

## Total sediment transport in overland flow structs and functions
@get_units @with_kw struct SedimentConcentrationsRiverVars{T}
    # Total sediment concentration in the river
    total::Vector{T} | "g m-3"
    # suspended sediemnt concentration in the river
    suspended::Vector{T} | "g m-3"
    # bed load sediment concentration in the river
    bed::Vector{T} | "g m-3"
end

function sediment_concentrations_river_vars(n)
    vars = SedimentConcentrationsRiverVars(;
        total = fill(mv, n),
        suspended = fill(mv, n),
        bed = fill(mv, n),
    )
    return vars
end

@get_units @with_kw struct SedimentConcentrationsRiverBC{T}
    # Discharge
    q::Vector{T} | "m3 s-1"
    waterlevel::Vector{T} | "m"
    # Clay load
    clay::Vector{T} | "g m-3"
    # Silt load
    silt::Vector{T} | "g m-3"
    # Sand load
    sand::Vector{T} | "g m-3"
    # Small aggregates load
    sagg::Vector{T} | "g m-3"
    # Large aggregates load
    lagg::Vector{T} | "g m-3"
    # Gravel load
    gravel::Vector{T} | "g m-3"
end

function sediment_concentrations_river_bc(n)
    bc = SedimentConcentrationsRiverBC(;
        q = fill(mv, n),
        waterlevel = fill(mv, n),
        clay = fill(mv, n),
        silt = fill(mv, n),
        sand = fill(mv, n),
        sagg = fill(mv, n),
        lagg = fill(mv, n),
        gravel = fill(mv, n),
    )
    return bc
end

# Common parameters for transport capacity models
@get_units @with_kw struct SedimentConcentrationsRiverParameters{T}
    # Clay mean diameter
    dm_clay::Vector{T} | "µm"
    # Silt mean diameter
    dm_silt::Vector{T} | "µm"
    # Sand mean diameter
    dm_sand::Vector{T} | "µm"
    # Small aggregates mean diameter
    dm_sagg::Vector{T} | "µm"
    # Large aggregates mean diameter
    dm_lagg::Vector{T} | "µm"
    # Gravel mean diameter
    dm_gravel::Vector{T} | "µm"
end

function initialize_sediment_concentrations_river_params(nc, config, inds)
    dm_clay = ncread(
        nc,
        config,
        "lateral.river.concentrations.parameters.dm_clay";
        sel = inds,
        defaults = 2.0,
        type = Float,
    )
    dm_silt = ncread(
        nc,
        config,
        "lateral.river.concentrations.parameters.dm_silt";
        sel = inds,
        defaults = 10.0,
        type = Float,
    )
    dm_sand = ncread(
        nc,
        config,
        "lateral.river.concentrations.parameters.dm_sand";
        sel = inds,
        defaults = 200.0,
        type = Float,
    )
    dm_sagg = ncread(
        nc,
        config,
        "lateral.river.concentrations.parameters.dm_sagg";
        sel = inds,
        defaults = 30.0,
        type = Float,
    )
    dm_lagg = ncread(
        nc,
        config,
        "lateral.river.concentrations.parameters.dm_lagg";
        sel = inds,
        defaults = 500.0,
        type = Float,
    )
    dm_gravel = ncread(
        nc,
        config,
        "lateral.river.concentrations.parameters.dm_gravel";
        sel = inds,
        defaults = 2000.0,
        type = Float,
    )
    conc_parameters = SedimentConcentrationsRiverParameters(;
        dm_clay = dm_clay,
        dm_silt = dm_silt,
        dm_sand = dm_sand,
        dm_sagg = dm_sagg,
        dm_lagg = dm_lagg,
        dm_gravel = dm_gravel,
    )

    return conc_parameters
end

@get_units @with_kw struct SedimentConcentrationsRiverModel{T} <:
                           AbstractSedimentConcentrationsRiverModel
    boundary_conditions::SedimentConcentrationsRiverBC{T} | "-"
    parameters::SedimentConcentrationsRiverParameters{T} | "-"
    variables::SedimentConcentrationsRiverVars{T} | "-"
end

function initialize_sediment_concentrations_river_model(nc, config, inds)
    n = length(inds)
    vars = sediment_concentrations_river_vars(n)
    params = initialize_sediment_concentrations_river_params(nc, config, inds)
    bc = sediment_concentrations_river_bc(n)
    model = SedimentConcentrationsRiverModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::SedimentConcentrationsRiverModel, geometry::RiverGeometry, ts)
    (; q, waterlevel, clay, silt, sand, sagg, lagg, gravel) = model.boundary_conditions
    (; dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg, dm_gravel) = model.parameters
    (; total, suspended, bed) = model.variables

    zeros = zeros(Float, length(q))
    # Conversion from load [ton] to concentration for rivers [mg/L]
    toconc .= ifelse.(q .> 0.0, 1e6 ./ (q .* ts), zeros)

    # Differentiation of bed and suspended load using Rouse number for suspension
    # threshold diameter between bed load and mixed load using Rouse number
    dbedf .=
        1e3 .*
        (2.5 .* 3600 .* 0.41 ./ 411 .* (9.81 .* waterlevel .* geometry.slope) .^ 0.5) .^ 0.5
    # threshold diameter between suspended load and mixed load using Rouse number
    dsuspf .=
        1e3 .*
        (1.2 .* 3600 .* 0.41 ./ 411 .* (9.81 .* waterlevel .* geometry.slope) .^ 0.5) .^ 0.5

    # Rouse with diameter
    SSclay .=
        ifelse.(dm_clay .<= dsuspf, clay, ifelse.(dm_clay .<= dbedf, clay ./ 2, zeros))
    SSsilt .=
        ifelse.(dm_silt .<= dsuspf, silt, ifelse.(dm_silt .<= dbedf, silt ./ 2, zeros))
    SSsand .=
        ifelse.(dm_sand .<= dsuspf, sand, ifelse.(dm_sand .<= dbedf, sand ./ 2, zeros))
    SSsagg .=
        ifelse.(dm_sagg .<= dsuspf, sagg, ifelse.(dm_sagg .<= dbedf, sagg ./ 2, zeros))
    SSlagg .=
        ifelse.(dm_lagg .<= dsuspf, lagg, ifelse.(dm_lagg .<= dbedf, lagg ./ 2, zeros))
    SSgrav .=
        ifelse.(
            dm_gravel .<= dsuspf,
            gravel,
            ifelse.(dm_gravel .<= dbedf, gravel ./ 2, zeros),
        )

    SS .= SSclay .+ SSsilt .+ SSsagg .+ SSsand .+ SSlagg .+ SSgrav
    Bedload .= (clay .+ silt .+ sagg .+ sand .+ lagg .+ gravel) .- SS

    suspended .= SS .* toconc
    bed .= Bedload .* toconc
    total .= (clay .+ silt .+ sagg .+ sand .+ lagg .+ gravel) .* toconc
end