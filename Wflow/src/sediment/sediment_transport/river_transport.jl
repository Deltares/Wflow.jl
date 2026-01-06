abstract type AbstractSedimentRiverTransportModel end

"Struct to store river sediment transport model variables"
@with_kw struct SedimentRiverTransportVariables
    n::Int
    # Sediment flux [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    clay_rate::Vector{Float64} = zeros(n)
    silt_rate::Vector{Float64} = zeros(n)
    sand_rate::Vector{Float64} = zeros(n)
    sagg_rate::Vector{Float64} = zeros(n)
    lagg_rate::Vector{Float64} = zeros(n)
    gravel_rate::Vector{Float64} = zeros(n)
    # Total Sediment deposition rate [t dt-1]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total sediment erosion rate (from store + direct river bed/bank) [t dt-1]
    erosion::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sediment / particle left in the cell [t] - states
    leftover_clay::Vector{Float64} = zeros(n)
    leftover_silt::Vector{Float64} = zeros(n)
    leftover_sand::Vector{Float64} = zeros(n)
    leftover_sagg::Vector{Float64} = zeros(n)
    leftover_lagg::Vector{Float64} = zeros(n)
    leftover_gravel::Vector{Float64} = zeros(n)
    # Sediment / particle stored on the river bed after deposition [t] -states
    store_clay::Vector{Float64} = zeros(n)
    store_silt::Vector{Float64} = zeros(n)
    store_sand::Vector{Float64} = zeros(n)
    store_sagg::Vector{Float64} = zeros(n)
    store_lagg::Vector{Float64} = zeros(n)
    store_gravel::Vector{Float64} = zeros(n)
end

"Struct to store river sediment transport model boundary conditions"
@with_kw struct SedimentRiverTransportBC
    n::Int
    # Waterlevel [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n)
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity of the flow [t dt-1]
    transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sediment input rate from land erosion [t dt-1]
    erosion_land_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    erosion_land_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    erosion_land_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    erosion_land_sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    erosion_land_lagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sediment erosion rate from direct river erosion [t dt-1]
    potential_erosion_river_bed::Vector{Float64} = fill(MISSING_VALUE, n)
    potential_erosion_river_bank::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store river sediment transport model parameters"
@with_kw struct SedimentRiverTransportParameters
    # River bed/bank content clay [-]
    clay_fraction::Vector{Float64}
    # River bed/bank content silt [-]
    silt_fraction::Vector{Float64}
    # River bed/bank content sand [-]
    sand_fraction::Vector{Float64}
    # River bed/bank content gravel [-]
    gravel_fraction::Vector{Float64}
    # Clay mean diameter [µm]
    dm_clay::Vector{Float64}
    # Silt mean diameter [µm]
    dm_silt::Vector{Float64}
    # Sand mean diameter [µm]
    dm_sand::Vector{Float64}
    # Small aggregates mean diameter [µm]
    dm_sagg::Vector{Float64}
    # Large aggregates mean diameter [µm]
    dm_lagg::Vector{Float64}
    # Gravel mean diameter [µm]
    dm_gravel::Vector{Float64}
    # Reservoir outlets [-]
    reservoir_outlet::Vector{Bool}
    # Reservoir area [m²]
    reservoir_area::Vector{Float64}
    # Reservoir trapping efficiency [-]
    reservoir_trapping_efficiency::Vector{Float64}
end

"Initialize river sediment transport model parameters"
function SedimentRiverTransportParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    clay_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_clay__mass_fraction";
        sel = indices,
        defaults = 0.15,
        type = Float64,
    )
    silt_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_silt__mass_fraction";
        sel = indices,
        defaults = 0.65,
        type = Float64,
    )
    sand_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_sand__mass_fraction";
        sel = indices,
        defaults = 0.15,
        type = Float64,
    )
    gravel_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_gravel__mass_fraction";
        sel = indices,
        defaults = 0.05,
        type = Float64,
    )
    # Check that river fractions sum to 1
    river_fractions = clay_fraction + silt_fraction + sand_fraction + gravel_fraction
    if any(abs.(river_fractions .- 1.0) .> 1e-3)
        error("Particle fractions in the river bed must sum to 1")
    end
    dm_clay = ncread(
        dataset,
        config,
        "clay__mean_diameter";
        sel = indices,
        defaults = 2.0,
        type = Float64,
    )
    dm_silt = ncread(
        dataset,
        config,
        "silt__mean_diameter";
        sel = indices,
        defaults = 10.0,
        type = Float64,
    )
    dm_sand = ncread(
        dataset,
        config,
        "sand__mean_diameter";
        sel = indices,
        defaults = 200.0,
        type = Float64,
    )
    dm_sagg = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter";
        sel = indices,
        defaults = 30.0,
        type = Float64,
    )
    dm_lagg = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter";
        sel = indices,
        defaults = 500.0,
        type = Float64,
    )
    dm_gravel = ncread(
        dataset,
        config,
        "gravel__mean_diameter";
        sel = indices,
        defaults = 2000.0,
        type = Float64,
    )

    # Reservoirs
    if config.model.reservoir__flag
        reservoir_outlet = ncread(
            dataset,
            config,
            "reservoir_location__count";
            sel = indices,
            type = Float64,
            fill = 0,
        )
        reservoir_area = ncread(
            dataset,
            config,
            "reservoir_surface__area";
            optional = false,
            sel = indices,
            type = Float64,
            fill = 0.0,
        )
        reservoir_trapping_efficiency = ncread(
            dataset,
            config,
            "reservoir_water_sediment__bedload_trapping_efficiency";
            sel = indices,
            type = Float64,
            defaults = 1.0,
            fill = 0.0,
        )
    else
        reservoir_outlet = zeros(n)
        reservoir_area = zeros(n)
        reservoir_trapping_efficiency = zeros(n)
    end

    river_parameters = SedimentRiverTransportParameters(;
        clay_fraction,
        silt_fraction,
        sand_fraction,
        gravel_fraction,
        dm_clay,
        dm_silt,
        dm_sand,
        dm_sagg,
        dm_lagg,
        dm_gravel,
        reservoir_outlet = reservoir_outlet .> 0,
        reservoir_area,
        reservoir_trapping_efficiency,
    )

    return river_parameters
end

"Struct to store river sediment transport model"
@with_kw struct SedimentRiverTransportModel <: AbstractSedimentRiverTransportModel
    n::Int
    boundary_conditions::SedimentRiverTransportBC = SedimentRiverTransportBC(; n)
    parameters::SedimentRiverTransportParameters
    variables::SedimentRiverTransportVariables = SedimentRiverTransportVariables(; n)
end

"Initialize river sediment transport model"
function SedimentRiverTransportModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)

    parameters = SedimentRiverTransportParameters(dataset, config, indices)
    model = SedimentRiverTransportModel(; n, parameters)
    return model
end

"Update boundary conditions for river sediment transport model"
function update_boundary_conditions!(
    model::SedimentRiverTransportModel,
    hydrological_forcing::HydrologicalForcing,
    transport_capacity_model::AbstractTransportCapacityModel,
    to_river_model::SedimentToRiverDifferentiationModel,
    potential_erosion_model::AbstractRiverErosionModel,
    indices_riv::Vector{Int},
)
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

    # Hydrological forcing
    (; q_river, waterlevel_river) = hydrological_forcing
    @. q = q_river
    @. waterlevel = waterlevel_river
    # Transport capacity
    @. transport_capacity = transport_capacity_model.variables.sediment_transport_capacity
    # Input from soil erosion
    (; clay_rate, silt_rate, sand_rate, sagg_rate, lagg_rate) = to_river_model.variables
    map!(i -> clay_rate[i], erosion_land_clay, indices_riv)
    map!(i -> silt_rate[i], erosion_land_silt, indices_riv)
    map!(i -> sand_rate[i], erosion_land_sand, indices_riv)
    map!(i -> sagg_rate[i], erosion_land_sagg, indices_riv)
    map!(i -> lagg_rate[i], erosion_land_lagg, indices_riv)
    # Maximum direct river bed/bank erosion
    @. potential_erosion_river_bed = potential_erosion_model.variables.bed
    @. potential_erosion_river_bank = potential_erosion_model.variables.bank
end

"Update river sediment transport model for a single timestep"
function update!(model::SedimentRiverTransportModel, domain::DomainRiver, dt::Float64)
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
    (;
        clay_fraction,
        silt_fraction,
        sand_fraction,
        gravel_fraction,
        dm_clay,
        dm_silt,
        dm_sand,
        dm_sagg,
        dm_lagg,
        dm_gravel,
        reservoir_outlet,
        reservoir_area,
        reservoir_trapping_efficiency,
    ) = model.parameters
    (;
        sediment_rate,
        clay_rate,
        silt_rate,
        sand_rate,
        sagg_rate,
        lagg_rate,
        gravel_rate,
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

    (; graph, order) = domain.network
    (; slope, flow_width, flow_length, reservoir_coverage) = domain.parameters

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
                if clay_rate[i] >= 0.0 # avoid NaN from upstream non-river cells
                    input_clay += clay_rate[i]
                    input_silt += silt_rate[i]
                    input_sand += sand_rate[i]
                    input_sagg += sagg_rate[i]
                    input_lagg += lagg_rate[i]
                    input_gravel += gravel_rate[i]
                end
            end
        end

        input_sediment =
            input_clay + input_silt + input_sand + input_sagg + input_lagg + input_gravel

        ### River erosion ###
        # Erosion only if the load is below the transport capacity of the flow.
        sediment_need = max(transport_capacity[v] - input_sediment, 0.0)
        # No erosion in reservoirs
        if reservoir_coverage[v]
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
            if (potential_erosion_river_bank[v] + potential_erosion_river_bed[v] > 0.0)
                RTEbank =
                    potential_erosion_river_bank[v] /
                    (potential_erosion_river_bank[v] + potential_erosion_river_bed[v])
            else
                RTEbank = 0.0
            end
            RTEbed = 1.0 - RTEbank

            # Actual bed and bank erosion
            erosion_bank = min(RTEbank * effsediment_need, potential_erosion_river_bank[v])
            erosion_bed = min(RTEbed * effsediment_need, potential_erosion_river_bed[v])
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
        # Different deposition if reservoir outlet or river
        if reservoir_outlet[v]
            # Deposition in reservoir outlets
            deposition_clay = reservoir_deposition_camp(
                (input_clay + erosion_clay),
                q[v],
                waterlevel[v],
                reservoir_area[v],
                reservoir_trapping_efficiency[v],
                dm_clay[v],
                slope[v],
            )
            deposition_silt = reservoir_deposition_camp(
                (input_silt + erosion_silt),
                q[v],
                waterlevel[v],
                reservoir_area[v],
                reservoir_trapping_efficiency[v],
                dm_silt[v],
                slope[v],
            )
            deposition_sand = reservoir_deposition_camp(
                (input_sand + erosion_sand),
                q[v],
                waterlevel[v],
                reservoir_area[v],
                reservoir_trapping_efficiency[v],
                dm_sand[v],
                slope[v],
            )
            deposition_sagg = reservoir_deposition_camp(
                (input_sagg + erosion_sagg),
                q[v],
                waterlevel[v],
                reservoir_area[v],
                reservoir_trapping_efficiency[v],
                dm_sagg[v],
                slope[v],
            )
            deposition_lagg = reservoir_deposition_camp(
                (input_lagg + erosion_lagg),
                q[v],
                waterlevel[v],
                reservoir_area[v],
                reservoir_trapping_efficiency[v],
                dm_lagg[v],
                slope[v],
            )
            deposition_gravel = reservoir_deposition_camp(
                (input_gravel + erosion_gravel),
                q[v],
                waterlevel[v],
                reservoir_area[v],
                reservoir_trapping_efficiency[v],
                dm_gravel[v],
                slope[v],
            )
        elseif reservoir_coverage[v]
            # No deposition in reservoir coverage, only at the outlets
            deposition_clay = 0.0
            deposition_silt = 0.0
            deposition_sand = 0.0
            deposition_sagg = 0.0
            deposition_lagg = 0.0
            deposition_gravel = 0.0
        else
            # Deposition in the river
            # From transport capacity exceedance
            excess_sediment = max(input_sediment - transport_capacity[v], 0.0)
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
                xs =
                    ifelse(q[v] > 0.0, 1.055 * flow_length[v] / (q[v] / flow_width[v]), 0.0)
                xclay =
                    min(1.0, 1.0 - 1.0 / exp(xs * (STOKES_FACTOR * (dm_clay[v] / 1000)^2)))
                deposition_clay = xclay * (input_clay + erosion_clay)
                xsilt =
                    min(1.0, 1.0 - 1.0 / exp(xs * (STOKES_FACTOR * (dm_silt[v] / 1000)^2)))
                deposition_silt = xsilt * (input_silt + erosion_silt)
                xsand =
                    min(1.0, 1.0 - 1.0 / exp(xs * (STOKES_FACTOR * (dm_sand[v] / 1000)^2)))
                deposition_sand = xsand * (input_sand + erosion_sand)
                xsagg =
                    min(1.0, 1.0 - 1.0 / exp(xs * (STOKES_FACTOR * (dm_sagg[v] / 1000)^2)))
                deposition_sagg = xsagg * (input_sagg + erosion_sagg)
                xlagg =
                    min(1.0, 1.0 - 1.0 / exp(xs * (STOKES_FACTOR * (dm_lagg[v] / 1000)^2)))
                deposition_lagg = xlagg * (input_lagg + erosion_lagg)
                xgrav = min(
                    1.0,
                    1.0 - 1.0 / exp(xs * (STOKES_FACTOR * (dm_gravel[v] / 1000)^2)),
                )
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
        if waterlevel[v] > 0.0
            fwaterout =
                min(q[v] * dt / (waterlevel[v] * flow_width[v] * flow_length[v]), 1.0)
        else
            fwaterout = 1.0
        end
        clay_rate[v] = fwaterout * (input_clay + erosion_clay - deposition_clay)
        silt_rate[v] = fwaterout * (input_silt + erosion_silt - deposition_silt)
        sand_rate[v] = fwaterout * (input_sand + erosion_sand - deposition_sand)
        sagg_rate[v] = fwaterout * (input_sagg + erosion_sagg - deposition_sagg)
        lagg_rate[v] = fwaterout * (input_lagg + erosion_lagg - deposition_lagg)
        gravel_rate[v] = fwaterout * (input_gravel + erosion_gravel - deposition_gravel)

        sediment_rate[v] =
            clay_rate[v] +
            silt_rate[v] +
            sand_rate[v] +
            sagg_rate[v] +
            lagg_rate[v] +
            gravel_rate[v]

        ### Leftover / mass balance ###
        # Sediment left in the cell [ton]
        leftover_clay[v] = input_clay + erosion_clay - deposition_clay - clay_rate[v]
        leftover_silt[v] = input_silt + erosion_silt - deposition_silt - silt_rate[v]
        leftover_sand[v] = input_sand + erosion_sand - deposition_sand - sand_rate[v]
        leftover_sagg[v] = input_sagg + erosion_sagg - deposition_sagg - sagg_rate[v]
        leftover_lagg[v] = input_lagg + erosion_lagg - deposition_lagg - lagg_rate[v]
        leftover_gravel[v] =
            input_gravel + erosion_gravel - deposition_gravel - gravel_rate[v]
    end
end

abstract type AbstractSedimentConcentrationsRiverModel end

"Struct to store river sediment concentrations model variables"
@with_kw struct SedimentConcentrationsRiverVariables
    n::Int
    # Total sediment concentration in the river [g m-3]
    total::Vector{Float64} = fill(MISSING_VALUE, n)
    # suspended sediemnt concentration in the river [g m-3]
    suspended::Vector{Float64} = fill(MISSING_VALUE, n)
    # bed load sediment concentration in the river [g m-3]
    bed::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store river sediment concentrations model boundary conditions"
@with_kw struct SedimentConcentrationsRiverBC
    n::Int
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n)
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n) # [m]
    # Clay load [g m-3]
    clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt load [g m-3]
    silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand load [g m-3]
    sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates load [g m-3]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates load [g m-3]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Gravel load [g m-3]
    gravel::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store river sediment concentrations model parameters"
@with_kw struct SedimentConcentrationsRiverParameters
    # Clay mean diameter [µm]
    dm_clay::Vector{Float64}
    # Silt mean diameter [µm]
    dm_silt::Vector{Float64}
    # Sand mean diameter [µm]
    dm_sand::Vector{Float64}
    # Small aggregates mean diameter [µm]
    dm_sagg::Vector{Float64}
    # Large aggregates mean diameter [µm]
    dm_lagg::Vector{Float64}
    # Gravel mean diameter [µm]
    dm_gravel::Vector{Float64}
end

"Initialize river sediment concentrations model parameters"
function SedimentConcentrationsRiverParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    dm_clay = ncread(
        dataset,
        config,
        "clay__mean_diameter";
        sel = indices,
        defaults = 2.0,
        type = Float64,
    )
    dm_silt = ncread(
        dataset,
        config,
        "silt__mean_diameter";
        sel = indices,
        defaults = 10.0,
        type = Float64,
    )
    dm_sand = ncread(
        dataset,
        config,
        "sand__mean_diameter";
        sel = indices,
        defaults = 200.0,
        type = Float64,
    )
    dm_sagg = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter";
        sel = indices,
        defaults = 30.0,
        type = Float64,
    )
    dm_lagg = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter";
        sel = indices,
        defaults = 500.0,
        type = Float64,
    )
    dm_gravel = ncread(
        dataset,
        config,
        "gravel__mean_diameter";
        sel = indices,
        defaults = 2000.0,
        type = Float64,
    )
    conc_parameters = SedimentConcentrationsRiverParameters(;
        dm_clay,
        dm_silt,
        dm_sand,
        dm_sagg,
        dm_lagg,
        dm_gravel,
    )

    return conc_parameters
end

"Struct to store river sediment concentrations model"
@with_kw struct SedimentConcentrationsRiverModel <: AbstractSedimentConcentrationsRiverModel
    n::Int
    boundary_conditions::SedimentConcentrationsRiverBC = SedimentConcentrationsRiverBC(; n)
    parameters::SedimentConcentrationsRiverParameters
    variables::SedimentConcentrationsRiverVariables =
        SedimentConcentrationsRiverVariables(; n)
end

"Initialize river sediment concentrations model"
function SedimentConcentrationsRiverModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = SedimentConcentrationsRiverParameters(dataset, config, indices)
    model = SedimentConcentrationsRiverModel(; n, parameters)
    return model
end

"Update boundary conditions for river sediment concentrations model"
function update_boundary_conditions!(
    model::SedimentConcentrationsRiverModel,
    hydrological_forcing::HydrologicalForcing,
    sediment_flux_model::AbstractSedimentRiverTransportModel,
)
    (; q, waterlevel, clay, silt, sand, sagg, lagg, gravel) = model.boundary_conditions
    # Hydrological forcing
    (; q_river, waterlevel_river) = hydrological_forcing
    @. q = q_river
    @. waterlevel = waterlevel_river
    # Sediment flux per particle
    @. clay = sediment_flux_model.variables.clay_rate
    @. silt = sediment_flux_model.variables.silt_rate
    @. sand = sediment_flux_model.variables.sand_rate
    @. sagg = sediment_flux_model.variables.sagg_rate
    @. lagg = sediment_flux_model.variables.lagg_rate
    @. gravel = sediment_flux_model.variables.gravel_rate
end

function suspended_solid(dm, dsuspf, dbedf, substance)
    return if dm <= dsuspf
        substance
    elseif dm <= dbedf
        substance / 2
    else
        0.0
    end
end

"Update river sediment concentrations model for a single timestep"
function update!(
    model::SedimentConcentrationsRiverModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; q, waterlevel, clay, silt, sand, sagg, lagg, gravel) = model.boundary_conditions
    (; dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg, dm_gravel) = model.parameters
    (; total, suspended, bed) = model.variables
    (; slope) = parameters

    for (i, flow) in enumerate(q)
        if flow > 0
            # Differentiation of bed and suspended load using Rouse number for suspension
            # threshold diameter between bed load and mixed load using Rouse number
            common_term =
                0.41 * sqrt(GRAVITATIONAL_ACCELERATION * waterlevel[i] * slope[i]) /
                STOKES_FACTOR
            dbedf = 1e3 * sqrt(2.5 * common_term)
            # # threshold diameter between suspended load and mixed load using Rouse number
            dsuspf = 1e3 * sqrt(1.2 * common_term)

            # Rouse with diameter
            SSclay = suspended_solid(dm_clay[i], dsuspf, dbedf, clay[i])
            SSsilt = suspended_solid(dm_silt[i], dsuspf, dbedf, silt[i])
            SSsand = suspended_solid(dm_sand[i], dsuspf, dbedf, sand[i])
            SSsagg = suspended_solid(dm_sagg[i], dsuspf, dbedf, sagg[i])
            SSlagg = suspended_solid(dm_lagg[i], dsuspf, dbedf, lagg[i])
            SSgrav = suspended_solid(dm_gravel[i], dsuspf, dbedf, gravel[i])

            to_conc = 1e6 / (q[i] * dt)
            total_ = clay[i] + silt[i] + sagg[i] + sand[i] + lagg[i] + gravel[i]
            total[i] = total_ * to_conc

            SS = SSclay + SSsilt + SSsand + SSsagg + SSlagg + SSgrav
            suspended[i] = SS * to_conc
            bed[i] = (total_ - SS) * to_conc
        else
            suspended[i] = 0.0
            bed[i] = 0.0
            total[i] = 0.0
        end
    end
end
