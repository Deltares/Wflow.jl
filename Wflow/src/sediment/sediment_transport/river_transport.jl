abstract type AbstractSedimentRiverTransportModel end

"Struct to store river sediment transport model variables"
@with_kw struct SedimentRiverTransportVariables
    n_cells::Int
    # Sediment flux [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    clay_rate::Vector{Float64} = zeros(n_cells)
    silt_rate::Vector{Float64} = zeros(n_cells)
    sand_rate::Vector{Float64} = zeros(n_cells)
    sagg_rate::Vector{Float64} = zeros(n_cells)
    lagg_rate::Vector{Float64} = zeros(n_cells)
    gravel_rate::Vector{Float64} = zeros(n_cells)
    # Total Sediment deposition rate [t dt-1]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total sediment erosion rate (from store + direct river bed/bank) [t dt-1]
    erosion::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Sediment / particle left in the cell [t] - states
    leftover_clay::Vector{Float64} = zeros(n_cells)
    leftover_silt::Vector{Float64} = zeros(n_cells)
    leftover_sand::Vector{Float64} = zeros(n_cells)
    leftover_sagg::Vector{Float64} = zeros(n_cells)
    leftover_lagg::Vector{Float64} = zeros(n_cells)
    leftover_gravel::Vector{Float64} = zeros(n_cells)
    # Sediment / particle stored on the river bed after deposition [t] -states
    store_clay::Vector{Float64} = zeros(n_cells)
    store_silt::Vector{Float64} = zeros(n_cells)
    store_sand::Vector{Float64} = zeros(n_cells)
    store_sagg::Vector{Float64} = zeros(n_cells)
    store_lagg::Vector{Float64} = zeros(n_cells)
    store_gravel::Vector{Float64} = zeros(n_cells)
end

"Struct to store river sediment transport model boundary conditions"
@with_kw struct SedimentRiverTransportBC
    n_cells::Int
    # Waterlevel [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Transport capacity of the flow [t dt-1]
    transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Sediment input rate from land erosion [t dt-1]
    erosion_land_clay::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    erosion_land_silt::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    erosion_land_sand::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    erosion_land_sagg::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    erosion_land_lagg::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Sediment erosion rate from direct river erosion [t dt-1]
    potential_erosion_river_bed::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    potential_erosion_river_bank::Vector{Float64} = fill(MISSING_VALUE, n_cells)
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
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_cells = length(river_indices_2d)
    clay_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_clay__mass_fraction",
        SoilLossModel;
        sel = river_indices_2d,
    )
    silt_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_silt__mass_fraction",
        SoilLossModel;
        sel = river_indices_2d,
    )
    sand_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_sand__mass_fraction",
        SoilLossModel;
        sel = river_indices_2d,
    )
    gravel_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_gravel__mass_fraction",
        SoilLossModel;
        sel = river_indices_2d,
    )
    # Check that river fractions sum to 1
    river_fractions = clay_fraction + silt_fraction + sand_fraction + gravel_fraction
    if any(abs.(river_fractions .- 1.0) .> 1e-3)
        error("Particle fractions in the river bed must sum to 1")
    end
    dm_clay = ncread(
        dataset,
        config,
        "clay__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_silt = ncread(
        dataset,
        config,
        "silt__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_sand = ncread(
        dataset,
        config,
        "sand__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_sagg = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_lagg = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_gravel = ncread(
        dataset,
        config,
        "gravel__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )

    # Reservoirs
    if config.model.reservoir__flag
        reservoir_outlet = ncread(
            dataset,
            config,
            "reservoir_location__count",
            SoilLossModel;
            sel = river_indices_2d,
        )
        reservoir_area = ncread(
            dataset,
            config,
            "reservoir_surface__area",
            SoilLossModel;
            sel = river_indices_2d,
        )
        reservoir_trapping_efficiency = ncread(
            dataset,
            config,
            "reservoir_water_sediment__bedload_trapping_efficiency",
            SoilLossModel;
            sel = river_indices_2d,
        )
    else
        reservoir_outlet = zeros(n_cells)
        reservoir_area = zeros(n_cells)
        reservoir_trapping_efficiency = zeros(n_cells)
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
    n_cells::Int
    boundary_conditions::SedimentRiverTransportBC =
        SedimentRiverTransportBC(; n_cells)
    parameters::SedimentRiverTransportParameters
    variables::SedimentRiverTransportVariables =
        SedimentRiverTransportVariables(; n_cells)
end

"Initialize river sediment transport model"
function SedimentRiverTransportModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_cells = length(river_indices_2d)

    parameters = SedimentRiverTransportParameters(dataset, config, river_indices_2d)
    sediment_transport_model = SedimentRiverTransportModel(; n_cells, parameters)
    return sediment_transport_model
end

"Update boundary conditions for river sediment transport model"
function update_bc_river_sediment_transport_model!(
    sediment_transport_model::SedimentRiverTransportModel,
    hydrological_forcing::HydrologicalForcing,
    transport_capacity_model::AbstractTransportCapacityModel,
    sediment_to_river_model::SedimentToRiverDifferentiationModel,
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
    ) = sediment_transport_model.boundary_conditions

    # Hydrological forcing
    (; q_river, waterlevel_river) = hydrological_forcing
    @. q = q_river
    @. waterlevel = waterlevel_river
    # Transport capacity
    @. transport_capacity = transport_capacity_model.variables.sediment_transport_capacity
    # Input from soil erosion
    (; clay_rate, silt_rate, sand_rate, sagg_rate, lagg_rate) =
        sediment_to_river_model.variables
    map!(cell_idx -> clay_rate[cell_idx], erosion_land_clay, indices_riv)
    map!(cell_idx -> silt_rate[cell_idx], erosion_land_silt, indices_riv)
    map!(cell_idx -> sand_rate[cell_idx], erosion_land_sand, indices_riv)
    map!(cell_idx -> sagg_rate[cell_idx], erosion_land_sagg, indices_riv)
    map!(cell_idx -> lagg_rate[cell_idx], erosion_land_lagg, indices_riv)
    # Maximum direct river bed/bank erosion
    @. potential_erosion_river_bed = potential_erosion_model.variables.bed
    @. potential_erosion_river_bank = potential_erosion_model.variables.bank
end

"""
Calculate sediment input from leftover sediment, land erosion, and upstream contributions
"""
function compute_sediment_input(
    sediment_transport_model::SedimentRiverTransportModel,
    graph::DiGraph,
    cell_idx::Int,
)
    (; boundary_conditions, variables) = sediment_transport_model
    (;
        erosion_land_clay,
        erosion_land_silt,
        erosion_land_sand,
        erosion_land_sagg,
        erosion_land_lagg,
    ) = boundary_conditions
    (;
        clay_rate,
        silt_rate,
        sand_rate,
        sagg_rate,
        lagg_rate,
        gravel_rate,
        leftover_clay,
        leftover_silt,
        leftover_sand,
        leftover_sagg,
        leftover_lagg,
        leftover_gravel,
    ) = variables

    # Base input from leftover and land erosion
    input_clay = leftover_clay[cell_idx] + erosion_land_clay[cell_idx]
    input_silt = leftover_silt[cell_idx] + erosion_land_silt[cell_idx]
    input_sand = leftover_sand[cell_idx] + erosion_land_sand[cell_idx]
    input_sagg = leftover_sagg[cell_idx] + erosion_land_sagg[cell_idx]
    input_lagg = leftover_lagg[cell_idx] + erosion_land_lagg[cell_idx]
    input_gravel = leftover_gravel[cell_idx]

    # Add upstream contribution
    upstream_nodes = inneighbors(graph, cell_idx)
    for upstream_cell_idx in upstream_nodes
        if clay_rate[upstream_cell_idx] >= 0.0 # avoid NaN from upstream non-river cells
            input_clay += clay_rate[upstream_cell_idx]
            input_silt += silt_rate[upstream_cell_idx]
            input_sand += sand_rate[upstream_cell_idx]
            input_sagg += sagg_rate[upstream_cell_idx]
            input_lagg += lagg_rate[upstream_cell_idx]
            input_gravel += gravel_rate[upstream_cell_idx]
        end
    end

    return (input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel)
end

"""
Calculate direct river bed/bank erosion based on sediment need
"""
function compute_direct_river_erosion(
    sediment_transport_model::SedimentRiverTransportModel,
    sediment_need::Float64,
    store_sediment::Float64,
    cell_idx::Int,
)
    (; potential_erosion_river_bed, potential_erosion_river_bank) =
        sediment_transport_model.boundary_conditions
    (; clay_fraction, silt_fraction, sand_fraction, gravel_fraction) =
        sediment_transport_model.parameters

    if sediment_need <= store_sediment
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    # Effective sediment needed from river bed and bank erosion [ton]
    effsediment_need = sediment_need - store_sediment

    # Relative potential erosion rates of the bed and the bank [-]
    if (
        potential_erosion_river_bank[cell_idx] +
        potential_erosion_river_bed[cell_idx] > 0.0
    )
        RTEbank =
            potential_erosion_river_bank[cell_idx] / (
                potential_erosion_river_bank[cell_idx] +
                potential_erosion_river_bed[cell_idx]
            )
    else
        RTEbank = 0.0
    end
    RTEbed = 1.0 - RTEbank

    # Actual bed and bank erosion
    erosion_bank =
        min(RTEbank * effsediment_need, potential_erosion_river_bank[cell_idx])
    erosion_bed =
        min(RTEbed * effsediment_need, potential_erosion_river_bed[cell_idx])
    erosion_river = erosion_bank + erosion_bed

    # Per particle
    erosion_clay = erosion_river * clay_fraction[cell_idx]
    erosion_silt = erosion_river * silt_fraction[cell_idx]
    erosion_sand = erosion_river * sand_fraction[cell_idx]
    erosion_gravel = erosion_river * gravel_fraction[cell_idx]
    # No small and large aggregates in the river bed/bank
    erosion_sagg = 0.0
    erosion_lagg = 0.0

    return (
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_sagg,
        erosion_lagg,
        erosion_gravel,
    )
end

"""
Calculate erosion from previously deposited sediment store
"""
function compute_store_erosion!(
    variables::SedimentRiverTransportVariables,
    sediment_need::Float64,
    cell_idx::Int,
)
    (; store_clay, store_silt, store_sand, store_sagg, store_lagg, store_gravel) = variables

    erosion_clay, sediment_need =
        river_erosion_store!(store_clay, sediment_need, cell_idx)
    erosion_silt, sediment_need =
        river_erosion_store!(store_silt, sediment_need, cell_idx)
    erosion_sagg, sediment_need =
        river_erosion_store!(store_sagg, sediment_need, cell_idx)
    erosion_sand, sediment_need =
        river_erosion_store!(store_sand, sediment_need, cell_idx)
    erosion_lagg, sediment_need =
        river_erosion_store!(store_lagg, sediment_need, cell_idx)
    erosion_gravel, _ = river_erosion_store!(store_gravel, sediment_need, cell_idx)

    return (
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_sagg,
        erosion_lagg,
        erosion_gravel,
    )
end

"""
Calculate sediment deposition in reservoir outlets using Camp's formula
"""
function compute_reservoir_deposition(
    sediment_transport_model::SedimentRiverTransportModel,
    domain_parameters::RiverParameters,
    input_particles::NTuple{6, Float64},
    erosion_particles::NTuple{6, Float64},
    cell_idx::Int,
)
    (; boundary_conditions, parameters) = sediment_transport_model
    (; q, waterlevel) = boundary_conditions
    (;
        dm_clay,
        dm_silt,
        dm_sand,
        dm_sagg,
        dm_lagg,
        dm_gravel,
        reservoir_area,
        reservoir_trapping_efficiency,
    ) = parameters
    (; slope) = domain_parameters

    input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel =
        input_particles
    erosion_clay, erosion_silt, erosion_sand, erosion_sagg, erosion_lagg, erosion_gravel =
        erosion_particles

    get_deposition(input, erosion, dm) = reservoir_deposition_camp(
        input + erosion,
        q[cell_idx],
        waterlevel[cell_idx],
        reservoir_area[cell_idx],
        reservoir_trapping_efficiency[cell_idx],
        dm[cell_idx],
        slope[cell_idx],
    )

    deposition_clay = get_deposition(input_clay, erosion_clay, dm_clay)
    deposition_silt = get_deposition(input_silt, erosion_silt, dm_silt)
    deposition_sand = get_deposition(input_sand, erosion_sand, dm_sand)
    deposition_sagg = get_deposition(input_sagg, erosion_sagg, dm_sagg)
    deposition_lagg = get_deposition(input_lagg, erosion_lagg, dm_lagg)
    deposition_gravel = get_deposition(input_gravel, erosion_gravel, dm_gravel)

    return (
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
        deposition_gravel,
    )
end

function transport_capacity_deposition(excess_sediment, input, erosion)
    total = input + erosion
    if excess_sediment > total
        deposition = total
        excess_sediment -= total
    else
        deposition = excess_sediment
        excess_sediment = 0.0
    end
    return deposition, excess_sediment
end

"""
Calculate river deposition from transport capacity exceedance (gravel to clay priority)
"""
function compute_transport_capacity_deposition(
    excess_sediment::Float64,
    input_particles::NTuple{6, Float64},
    erosion_particles::NTuple{6, Float64},
)
    input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel =
        input_particles
    erosion_clay, erosion_silt, erosion_sand, erosion_sagg, erosion_lagg, erosion_gravel =
        erosion_particles

    # From largest to smallest particles
    deposition_gravel, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_gravel, erosion_gravel)
    deposition_lagg, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_lagg, erosion_lagg)
    deposition_sand, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_sand, erosion_sand)
    deposition_sagg, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_sagg, erosion_sagg)
    deposition_silt, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_silt, erosion_silt)
    deposition_clay, _ =
        transport_capacity_deposition(excess_sediment, input_clay, erosion_clay)

    return (
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
        deposition_gravel,
    )
end

function natural_deposition(xs::Float64, dm::Float64, input::Float64, erosion::Float64)
    x_material = min(1.0, 1.0 - exp(-xs * (STOKES_FACTOR * (dm / 1000)^2)))
    deposition = x_material * (input + erosion)
    return deposition
end

"""
Calculate natural river deposition using Einstein's formula (Stokes settling)
"""
function compute_natural_deposition(
    sediment_transport_model::SedimentRiverTransportModel,
    domain_parameters::RiverParameters,
    input_particles::NTuple{6, Float64},
    erosion_particles::NTuple{6, Float64},
    cell_idx::Int,
)
    (; boundary_conditions, parameters) = sediment_transport_model
    (; q) = boundary_conditions
    (; dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg, dm_gravel) = parameters
    (; flow_length, flow_width) = domain_parameters

    input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel =
        input_particles
    erosion_clay, erosion_silt, erosion_sand, erosion_sagg, erosion_lagg, erosion_gravel =
        erosion_particles

    # Particle fall velocity [m/s] from Stokes
    xs = ifelse(
        q[cell_idx] > 0.0,
        1.055 * flow_length[cell_idx] /
        (q[cell_idx] / flow_width[cell_idx]),
        0.0,
    )

    deposition_clay =
        natural_deposition(xs, dm_clay[cell_idx], input_clay, erosion_clay)
    deposition_silt =
        natural_deposition(xs, dm_silt[cell_idx], input_silt, erosion_silt)
    deposition_sand =
        natural_deposition(xs, dm_sand[cell_idx], input_sand, erosion_sand)
    deposition_sagg =
        natural_deposition(xs, dm_sagg[cell_idx], input_sagg, erosion_sagg)
    deposition_lagg =
        natural_deposition(xs, dm_lagg[cell_idx], input_lagg, erosion_lagg)
    deposition_gravel =
        natural_deposition(xs, dm_gravel[cell_idx], input_gravel, erosion_gravel)

    return (
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
        deposition_gravel,
    )
end

function water_outflow_fraction(waterlevel, q, flow_width, flow_length, dt)
    return if waterlevel > 0.0
        min(q * dt / (waterlevel * flow_width * flow_length), 1.0)
    else
        1.0
    end
end

function update_variables!(
    variables::SedimentRiverTransportVariables,
    input_particles::NTuple{6, Float64},
    erosion_particles::NTuple{6, Float64},
    deposition_particles::NTuple{6, Float64},
    fwaterout::Float64,
    cell_idx::Int,
)
    (;
        store_clay,
        store_silt,
        store_sand,
        store_sagg,
        store_lagg,
        store_gravel,
        deposition,
        erosion,
        sediment_rate,
        clay_rate,
        silt_rate,
        sand_rate,
        sagg_rate,
        lagg_rate,
        gravel_rate,
        leftover_clay,
        leftover_silt,
        leftover_sand,
        leftover_sagg,
        leftover_lagg,
        leftover_gravel,
    ) = variables

    (input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel) =
        input_particles
    (erosion_clay, erosion_silt, erosion_sand, erosion_sagg, erosion_lagg, erosion_gravel) =
        erosion_particles
    (
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
        deposition_gravel,
    ) = deposition_particles

    # Update the sediment store
    store_clay[cell_idx] += deposition_clay
    store_silt[cell_idx] += deposition_silt
    store_sand[cell_idx] += deposition_sand
    store_sagg[cell_idx] += deposition_sagg
    store_lagg[cell_idx] += deposition_lagg
    store_gravel[cell_idx] += deposition_gravel

    # Compute total erosion
    erosion[cell_idx] = sum(erosion_particles)

    # Compute total deposition
    deposition[cell_idx] = sum(deposition_particles)

    # Output loads
    clay_rate[cell_idx] = fwaterout * (input_clay + erosion_clay - deposition_clay)
    silt_rate[cell_idx] = fwaterout * (input_silt + erosion_silt - deposition_silt)
    sand_rate[cell_idx] = fwaterout * (input_sand + erosion_sand - deposition_sand)
    sagg_rate[cell_idx] = fwaterout * (input_sagg + erosion_sagg - deposition_sagg)
    lagg_rate[cell_idx] = fwaterout * (input_lagg + erosion_lagg - deposition_lagg)
    gravel_rate[cell_idx] =
        fwaterout * (input_gravel + erosion_gravel - deposition_gravel)

    sediment_rate[cell_idx] =
        clay_rate[cell_idx] +
        silt_rate[cell_idx] +
        sand_rate[cell_idx] +
        sagg_rate[cell_idx] +
        lagg_rate[cell_idx] +
        gravel_rate[cell_idx]

    # Leftover sediment for mass balance
    leftover_clay[cell_idx] =
        input_clay + erosion_clay - deposition_clay - clay_rate[cell_idx]
    leftover_silt[cell_idx] =
        input_silt + erosion_silt - deposition_silt - silt_rate[cell_idx]
    leftover_sand[cell_idx] =
        input_sand + erosion_sand - deposition_sand - sand_rate[cell_idx]
    leftover_sagg[cell_idx] =
        input_sagg + erosion_sagg - deposition_sagg - sagg_rate[cell_idx]
    leftover_lagg[cell_idx] =
        input_lagg + erosion_lagg - deposition_lagg - lagg_rate[cell_idx]
    leftover_gravel[cell_idx] =
        input_gravel + erosion_gravel - deposition_gravel - gravel_rate[cell_idx]

    return nothing
end

"Update river sediment transport model for a single timestep"
function update_sediment_river_transport_model!(
    sediment_transport_model::SedimentRiverTransportModel,
    domain::DomainRiver,
    dt::Float64,
)
    (; waterlevel, q, transport_capacity) = sediment_transport_model.boundary_conditions
    (; reservoir_outlet) = sediment_transport_model.parameters
    (; store_clay, store_silt, store_sand, store_sagg, store_lagg, store_gravel) =
        sediment_transport_model.variables

    (; graph, cell_order) = domain.network
    (; flow_width, flow_length, reservoir_coverage) = domain.parameters

    # Sediment transport - water balance in the river
    for cell_idx in cell_order
        ### Sediment input in the cell (left from previous timestep + from land + from upstream outflux) ###
        input_particles =
            compute_sediment_input(sediment_transport_model, graph, cell_idx)
        input_sediment = sum(input_particles)

        ### River erosion ###
        # Erosion only if the load is below the transport capacity of the flow.
        sediment_need = max(transport_capacity[cell_idx] - input_sediment, 0.0)
        # No erosion in reservoirs
        if reservoir_coverage[cell_idx]
            sediment_need = 0.0
        end

        # Available sediment stored from previous deposition
        store_sediment =
            store_clay[cell_idx] +
            store_silt[cell_idx] +
            store_sand[cell_idx] +
            store_sagg[cell_idx] +
            store_lagg[cell_idx] +
            store_gravel[cell_idx]

        # Direct erosion from the river bed/bank
        erosion_particles = compute_direct_river_erosion(
            sediment_transport_model,
            sediment_need,
            store_sediment,
            cell_idx,
        )

        # Erosion/degradation of the previously deposited sediment (from clay to gravel) [ton]
        store_erosion_particles = compute_store_erosion!(
            sediment_transport_model.variables,
            sediment_need,
            cell_idx,
        )

        # Update total erosion
        erosion_particles = erosion_particles .+ store_erosion_particles

        ### Deposition / settling ###

        # Different deposition if reservoir outlet or river
        deposition_particles = if reservoir_outlet[cell_idx]
            # Deposition in reservoir outlets
            compute_reservoir_deposition(
                sediment_transport_model,
                domain.parameters,
                input_particles,
                erosion_particles,
                cell_idx,
            )
        elseif reservoir_coverage[cell_idx]
            # No deposition in reservoir coverage, only at the outlets
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        else
            # Deposition in the river

            # From transport capacity exceedance
            excess_sediment = max(input_sediment - transport_capacity[cell_idx], 0.0)
            if excess_sediment > 0.0
                compute_transport_capacity_deposition(
                    excess_sediment,
                    input_particles,
                    erosion_particles,
                )
            else
                # Natural deposition from Einstein's formula (density controlled)
                compute_natural_deposition(
                    sediment_transport_model,
                    domain.parameters,
                    input_particles,
                    erosion_particles,
                    cell_idx,
                )
            end
        end

        fwaterout = water_outflow_fraction(
            waterlevel[cell_idx],
            q[cell_idx],
            flow_width[cell_idx],
            flow_length[cell_idx],
            dt,
        )

        update_variables!(
            sediment_transport_model.variables,
            input_particles,
            erosion_particles,
            deposition_particles,
            fwaterout,
            cell_idx,
        )
    end
    return nothing
end

abstract type AbstractSedimentConcentrationsRiverModel end

"Struct to store river sediment concentrations model variables"
@with_kw struct SedimentConcentrationsRiverVariables
    n_cells::Int
    # Total sediment concentration in the river [g m-3]
    total::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # suspended sediemnt concentration in the river [g m-3]
    suspended::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # bed load sediment concentration in the river [g m-3]
    bed::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct to store river sediment concentrations model boundary conditions"
@with_kw struct SedimentConcentrationsRiverBC
    n_cells::Int
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n_cells) # [m]
    # Clay load [g m-3]
    clay::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Silt load [g m-3]
    silt::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Sand load [g m-3]
    sand::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Small aggregates load [g m-3]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Large aggregates load [g m-3]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Gravel load [g m-3]
    gravel::Vector{Float64} = fill(MISSING_VALUE, n_cells)
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
    river_indices_2d::Vector{CartesianIndex{2}},
)
    dm_clay = ncread(
        dataset,
        config,
        "clay__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_silt = ncread(
        dataset,
        config,
        "silt__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_sand = ncread(
        dataset,
        config,
        "sand__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_sagg = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_lagg = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
    )
    dm_gravel = ncread(
        dataset,
        config,
        "gravel__mean_diameter",
        SoilLossModel;
        sel = river_indices_2d,
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
    n_cells::Int
    boundary_conditions::SedimentConcentrationsRiverBC =
        SedimentConcentrationsRiverBC(; n_cells)
    parameters::SedimentConcentrationsRiverParameters
    variables::SedimentConcentrationsRiverVariables =
        SedimentConcentrationsRiverVariables(; n_cells)
end

"Initialize river sediment concentrations model"
function SedimentConcentrationsRiverModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    n_cells = length(river_indices_2d)
    parameters = SedimentConcentrationsRiverParameters(dataset, config, river_indices_2d)
    sediment_concentrations_model =
        SedimentConcentrationsRiverModel(; n_cells, parameters)
    return sediment_concentrations_model
end

"Update boundary conditions for river sediment concentrations model"
function update_bc_river_sediment_concentration_model!(
    sediment_concentrations_model::SedimentConcentrationsRiverModel,
    hydrological_forcing::HydrologicalForcing,
    sediment_transport_model::AbstractSedimentRiverTransportModel,
)
    (; q, waterlevel, clay, silt, sand, sagg, lagg, gravel) =
        sediment_concentrations_model.boundary_conditions
    # Hydrological forcing
    (; q_river, waterlevel_river) = hydrological_forcing
    @. q = q_river
    @. waterlevel = waterlevel_river
    # Sediment flux per particle
    @. clay = sediment_transport_model.variables.clay_rate
    @. silt = sediment_transport_model.variables.silt_rate
    @. sand = sediment_transport_model.variables.sand_rate
    @. sagg = sediment_transport_model.variables.sagg_rate
    @. lagg = sediment_transport_model.variables.lagg_rate
    @. gravel = sediment_transport_model.variables.gravel_rate
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
function update_river_sediment_concentration_model!(
    sediment_transport_model::SedimentConcentrationsRiverModel,
    parameters::RiverParameters,
    dt::Float64,
)
    (; n_cells) = sediment_transport_model
    (; q, waterlevel, clay, silt, sand, sagg, lagg, gravel) =
        sediment_transport_model.boundary_conditions
    (; dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg, dm_gravel) =
        sediment_transport_model.parameters
    (; total, suspended, bed) = sediment_transport_model.variables
    (; slope) = parameters

    for cell_idx in 1:n_cells
        flow = q[cell_idx]
        if flow > 0
            # Differentiation of bed and suspended load using Rouse number for suspension
            # threshold diameter between bed load and mixed load using Rouse number
            common_term =
                0.41 * sqrt(
                    GRAVITATIONAL_ACCELERATION *
                    waterlevel[cell_idx] *
                    slope[cell_idx],
                ) / STOKES_FACTOR
            dbedf = 1e3 * sqrt(2.5 * common_term)
            # # threshold diameter between suspended load and mixed load using Rouse number
            dsuspf = 1e3 * sqrt(1.2 * common_term)

            # Rouse with diameter
            SSclay = suspended_solid(
                dm_clay[cell_idx],
                dsuspf,
                dbedf,
                clay[cell_idx],
            )
            SSsilt = suspended_solid(
                dm_silt[cell_idx],
                dsuspf,
                dbedf,
                silt[cell_idx],
            )
            SSsand = suspended_solid(
                dm_sand[cell_idx],
                dsuspf,
                dbedf,
                sand[cell_idx],
            )
            SSsagg = suspended_solid(
                dm_sagg[cell_idx],
                dsuspf,
                dbedf,
                sagg[cell_idx],
            )
            SSlagg = suspended_solid(
                dm_lagg[cell_idx],
                dsuspf,
                dbedf,
                lagg[cell_idx],
            )
            SSgrav = suspended_solid(
                dm_gravel[cell_idx],
                dsuspf,
                dbedf,
                gravel[cell_idx],
            )

            to_conc = 1e6 / (q[cell_idx] * dt)
            total_ =
                clay[cell_idx] +
                silt[cell_idx] +
                sagg[cell_idx] +
                sand[cell_idx] +
                lagg[cell_idx] +
                gravel[cell_idx]
            total[cell_idx] = total_ * to_conc

            SS = SSclay + SSsilt + SSsand + SSsagg + SSlagg + SSgrav
            suspended[cell_idx] = SS * to_conc
            bed[cell_idx] = (total_ - SS) * to_conc
        else
            suspended[cell_idx] = 0.0
            bed[cell_idx] = 0.0
            total[cell_idx] = 0.0
        end
    end
end
