abstract type AbstractSedimentRiverTransportModel end
abstract type AbstractSedimentConcentrationsRiverModel end

"Struct to store river sediment transport model variables"
@with_kw struct SedimentRiverTransportVariables
    n_river_cells::Int
    # Sediment flux [kg s⁻¹]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    clay_rate::Vector{Float64} = zeros(n_river_cells)
    silt_rate::Vector{Float64} = zeros(n_river_cells)
    sand_rate::Vector{Float64} = zeros(n_river_cells)
    small_aggregates_rate::Vector{Float64} = zeros(n_river_cells)
    large_aggregates_rate::Vector{Float64} = zeros(n_river_cells)
    gravel_rate::Vector{Float64} = zeros(n_river_cells)
    # Total Sediment deposition rate [kg s⁻¹]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Total sediment erosion rate (from store + direct river bed/bank) [kg s⁻¹]
    erosion::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Sediment / particle left in the cell [kg] - states
    leftover_clay::Vector{Float64} = zeros(n_river_cells)
    leftover_silt::Vector{Float64} = zeros(n_river_cells)
    leftover_sand::Vector{Float64} = zeros(n_river_cells)
    leftover_small_aggregates::Vector{Float64} = zeros(n_river_cells)
    leftover_large_aggregates::Vector{Float64} = zeros(n_river_cells)
    leftover_gravel::Vector{Float64} = zeros(n_river_cells)
    # Sediment / particle stored on the river bed after deposition [kg] -states
    store_clay::Vector{Float64} = zeros(n_river_cells)
    store_silt::Vector{Float64} = zeros(n_river_cells)
    store_sand::Vector{Float64} = zeros(n_river_cells)
    store_small_aggregates::Vector{Float64} = zeros(n_river_cells)
    store_large_aggregates::Vector{Float64} = zeros(n_river_cells)
    store_gravel::Vector{Float64} = zeros(n_river_cells)
end

"Struct to store river sediment transport model boundary conditions"
@with_kw struct SedimentRiverTransportBC
    n_river_cells::Int
    # Waterlevel [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Transport capacity of the flow [kg s⁻¹]
    transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Sediment input rate from land erosion [kg s⁻¹]
    erosion_land_clay::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    erosion_land_silt::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    erosion_land_sand::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    erosion_land_small_aggregates::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    erosion_land_large_aggregates::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Sediment erosion rate from direct river erosion [kg s⁻¹]
    potential_erosion_river_bed::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    potential_erosion_river_bank::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
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
    # Clay mean diameter [m]
    median_diameter_clay::Vector{Float64}
    # Silt mean diameter [m]
    median_diameter_silt::Vector{Float64}
    # Sand mean diameter [m]
    median_diameter_sand::Vector{Float64}
    # Small aggregates mean diameter [m]
    median_diameter_small_aggregates::Vector{Float64}
    # Large aggregates mean diameter [m]
    median_diameter_large_aggregates::Vector{Float64}
    # Gravel mean diameter [m]
    median_diameter_gravel::Vector{Float64}
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
    n_river_cells = length(indices)
    clay_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_clay__mass_fraction",
        SoilLossModel;
        sel = indices,
    )
    silt_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_silt__mass_fraction",
        SoilLossModel;
        sel = indices,
    )
    sand_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_sand__mass_fraction",
        SoilLossModel;
        sel = indices,
    )
    gravel_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_gravel__mass_fraction",
        SoilLossModel;
        sel = indices,
    )
    # Check that river fractions sum to 1
    river_fractions = clay_fraction + silt_fraction + sand_fraction + gravel_fraction
    if any(abs.(river_fractions .- 1.0) .> 1e-3)
        error("Particle fractions in the river bed must sum to 1")
    end
    median_diameter_clay =
        ncread(dataset, config, "clay__mean_diameter", SoilLossModel; sel = indices)
    median_diameter_silt =
        ncread(dataset, config, "silt__mean_diameter", SoilLossModel; sel = indices)
    median_diameter_sand =
        ncread(dataset, config, "sand__mean_diameter", SoilLossModel; sel = indices)
    median_diameter_small_aggregates = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter",
        SoilLossModel;
        sel = indices,
    )
    median_diameter_large_aggregates = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter",
        SoilLossModel;
        sel = indices,
    )
    median_diameter_gravel =
        ncread(dataset, config, "gravel__mean_diameter", SoilLossModel; sel = indices)

    # Reservoirs
    if config.model.reservoir__flag
        reservoir_outlet = ncread(
            dataset,
            config,
            "reservoir_location__count",
            SoilLossModel;
            sel = indices,
        )
        reservoir_area =
            ncread(dataset, config, "reservoir_surface__area", SoilLossModel; sel = indices)
        reservoir_trapping_efficiency = ncread(
            dataset,
            config,
            "reservoir_water_sediment__bedload_trapping_efficiency",
            SoilLossModel;
            sel = indices,
        )
    else
        reservoir_outlet = zeros(n_river_cells)
        reservoir_area = zeros(n_river_cells)
        reservoir_trapping_efficiency = zeros(n_river_cells)
    end

    river_parameters = SedimentRiverTransportParameters(;
        clay_fraction,
        silt_fraction,
        sand_fraction,
        gravel_fraction,
        median_diameter_clay,
        median_diameter_silt,
        median_diameter_sand,
        median_diameter_small_aggregates,
        median_diameter_large_aggregates,
        median_diameter_gravel,
        reservoir_outlet = reservoir_outlet .> 0,
        reservoir_area,
        reservoir_trapping_efficiency,
    )

    return river_parameters
end

"Struct to store river sediment transport model"
@with_kw struct SedimentRiverTransportModel <: AbstractSedimentRiverTransportModel
    n_river_cells::Int
    boundary_conditions::SedimentRiverTransportBC =
        SedimentRiverTransportBC(; n_river_cells)
    parameters::SedimentRiverTransportParameters
    variables::SedimentRiverTransportVariables =
        SedimentRiverTransportVariables(; n_river_cells)
end

"Initialize river sediment transport model"
function SedimentRiverTransportModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n_river_cells = length(indices)
    parameters = SedimentRiverTransportParameters(dataset, config, indices)
    sediment_transport_model = SedimentRiverTransportModel(; n_river_cells, parameters)
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
        erosion_land_small_aggregates,
        erosion_land_large_aggregates,
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
    (; clay_rate, silt_rate, sand_rate, small_aggregates_rate, large_aggregates_rate) =
        sediment_to_river_model.variables
    map!(i -> clay_rate[i], erosion_land_clay, indices_riv)
    map!(i -> silt_rate[i], erosion_land_silt, indices_riv)
    map!(i -> sand_rate[i], erosion_land_sand, indices_riv)
    map!(i -> small_aggregates_rate[i], erosion_land_small_aggregates, indices_riv)
    map!(i -> large_aggregates_rate[i], erosion_land_large_aggregates, indices_riv)
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
    dt::Float64,
    v::Int,
)
    (; boundary_conditions, variables) = sediment_transport_model
    (;
        erosion_land_clay,
        erosion_land_silt,
        erosion_land_sand,
        erosion_land_small_aggregates,
        erosion_land_large_aggregates,
    ) = boundary_conditions
    (;
        clay_rate,
        silt_rate,
        sand_rate,
        small_aggregates_rate,
        large_aggregates_rate,
        gravel_rate,
        leftover_clay,
        leftover_silt,
        leftover_sand,
        leftover_small_aggregates,
        leftover_large_aggregates,
        leftover_gravel,
    ) = variables

    # Base input from leftover and land erosion
    input_clay = leftover_clay[v] / dt + erosion_land_clay[v]
    input_silt = leftover_silt[v] / dt + erosion_land_silt[v]
    input_sand = leftover_sand[v] / dt + erosion_land_sand[v]
    input_sagg = leftover_small_aggregates[v] / dt + erosion_land_small_aggregates[v]
    input_lagg = leftover_large_aggregates[v] / dt + erosion_land_large_aggregates[v]
    input_gravel = leftover_gravel[v] / dt

    # Add upstream contribution
    upstream_nodes = inneighbors(graph, v)
    if !isempty(upstream_nodes)
        for river_idx in upstream_nodes
            if clay_rate[river_idx] >= 0.0 # avoid NaN from upstream non-river cells
                input_clay += clay_rate[river_idx]
                input_silt += silt_rate[river_idx]
                input_sand += sand_rate[river_idx]
                input_sagg += small_aggregates_rate[river_idx]
                input_lagg += large_aggregates_rate[river_idx]
                input_gravel += gravel_rate[river_idx]
            end
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
    dt::Float64,
    v::Int,
)
    (; potential_erosion_river_bed, potential_erosion_river_bank) =
        sediment_transport_model.boundary_conditions
    (; clay_fraction, silt_fraction, sand_fraction, gravel_fraction) =
        sediment_transport_model.parameters

    if sediment_need <= store_sediment / dt
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    # Effective sediment needed from river bed and bank erosion
    effsediment_need = sediment_need - store_sediment / dt

    # Relative potential erosion rates of the bed and the bank
    RTEbank = if (potential_erosion_river_bank[v] + potential_erosion_river_bed[v] > 0.0)
        potential_erosion_river_bank[v] /
        (potential_erosion_river_bank[v] + potential_erosion_river_bed[v])
    else
        0.0
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
    erosion_small_aggregates = 0.0
    erosion_large_aggregates = 0.0

    return (
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_small_aggregates,
        erosion_large_aggregates,
        erosion_gravel,
    )
end

"""
Calculate erosion from previously deposited sediment store
"""
function compute_store_erosion!(
    variables::SedimentRiverTransportVariables,
    sediment_need::Float64,
    dt::Float64,
    v::Int,
)
    (;
        store_clay,
        store_silt,
        store_sand,
        store_small_aggregates,
        store_large_aggregates,
        store_gravel,
    ) = variables

    erosion_clay, sediment_need = river_erosion_store!(store_clay, sediment_need, dt, v)
    erosion_silt, sediment_need = river_erosion_store!(store_silt, sediment_need, dt, v)
    erosion_small_aggregates, sediment_need =
        river_erosion_store!(store_small_aggregates, sediment_need, dt, v)
    erosion_sand, sediment_need = river_erosion_store!(store_sand, sediment_need, dt, v)
    erosion_large_aggregates, sediment_need =
        river_erosion_store!(store_large_aggregates, sediment_need, dt, v)
    erosion_gravel, _ = river_erosion_store!(store_gravel, sediment_need, dt, v)

    return (
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_small_aggregates,
        erosion_large_aggregates,
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
    v::Int,
)
    (; boundary_conditions, parameters) = sediment_transport_model
    (; q, waterlevel) = boundary_conditions
    (;
        median_diameter_clay,
        median_diameter_silt,
        median_diameter_sand,
        median_diameter_small_aggregates,
        median_diameter_large_aggregates,
        median_diameter_gravel,
        reservoir_area,
        reservoir_trapping_efficiency,
    ) = parameters
    (; slope) = domain_parameters

    input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel =
        input_particles
    erosion_clay,
    erosion_silt,
    erosion_sand,
    erosion_small_aggregates,
    erosion_large_aggregates,
    erosion_gravel = erosion_particles

    get_deposition(input, erosion, dm) = reservoir_deposition_camp(
        input + erosion,
        q[v],
        waterlevel[v],
        reservoir_area[v],
        reservoir_trapping_efficiency[v],
        dm[v],
        slope[v],
    )

    deposition_clay = get_deposition(input_clay, erosion_clay, median_diameter_clay)
    deposition_silt = get_deposition(input_silt, erosion_silt, median_diameter_silt)
    deposition_sand = get_deposition(input_sand, erosion_sand, median_diameter_sand)
    deposition_small_aggregates = get_deposition(
        input_sagg,
        erosion_small_aggregates,
        median_diameter_small_aggregates,
    )
    deposition_large_aggregates = get_deposition(
        input_lagg,
        erosion_large_aggregates,
        median_diameter_large_aggregates,
    )
    deposition_gravel = get_deposition(input_gravel, erosion_gravel, median_diameter_gravel)

    return (
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_small_aggregates,
        deposition_large_aggregates,
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

function compute_transport_capacity_deposition(
    excess_sediment,
    input_particles,
    erosion_particles,
)
    input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel =
        input_particles
    erosion_clay,
    erosion_silt,
    erosion_sand,
    erosion_small_aggregates,
    erosion_large_aggregates,
    erosion_gravel = erosion_particles

    # From largest to smallest particles
    deposition_gravel, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_gravel, erosion_gravel)
    deposition_large_aggregates, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_lagg, erosion_large_aggregates)
    deposition_sand, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_sand, erosion_sand)
    deposition_small_aggregates, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_sagg, erosion_small_aggregates)
    deposition_silt, excess_sediment =
        transport_capacity_deposition(excess_sediment, input_silt, erosion_silt)
    deposition_clay, _ =
        transport_capacity_deposition(excess_sediment, input_clay, erosion_clay)

    return (
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_small_aggregates,
        deposition_large_aggregates,
        deposition_gravel,
    )
end

function natural_deposition(xs::Float64, dm::Float64, input::Float64, erosion::Float64)
    x_material = min(1.0, 1.0 - exp(-xs * fall_velocity(dm)))
    return x_material * (input + erosion)
end

"""
Calculate natural river deposition using Einstein's formula (Stokes settling)
"""
function compute_natural_deposition(
    sediment_transport_model::SedimentRiverTransportModel,
    domain_parameters::RiverParameters,
    input_particles::NTuple{6, Float64},
    erosion_particles::NTuple{6, Float64},
    v::Int,
)
    (; boundary_conditions, parameters) = sediment_transport_model
    (; q) = boundary_conditions
    (;
        median_diameter_clay,
        median_diameter_silt,
        median_diameter_sand,
        median_diameter_small_aggregates,
        median_diameter_large_aggregates,
        median_diameter_gravel,
    ) = parameters
    (; flow_length, flow_width) = domain_parameters

    input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel =
        input_particles
    erosion_clay,
    erosion_silt,
    erosion_sand,
    erosion_small_aggregates,
    erosion_large_aggregates,
    erosion_gravel = erosion_particles

    # Particle fall velocity [m/s] from Stokes
    xs = ifelse(q[v] > 0.0, 1.055 * flow_length[v] / (q[v] / flow_width[v]), 0.0)

    deposition_clay =
        natural_deposition(xs, median_diameter_clay[v], input_clay, erosion_clay)
    deposition_silt =
        natural_deposition(xs, median_diameter_silt[v], input_silt, erosion_silt)
    deposition_sand =
        natural_deposition(xs, median_diameter_sand[v], input_sand, erosion_sand)
    deposition_small_aggregates = natural_deposition(
        xs,
        median_diameter_small_aggregates[v],
        input_sagg,
        erosion_small_aggregates,
    )
    deposition_large_aggregates = natural_deposition(
        xs,
        median_diameter_large_aggregates[v],
        input_lagg,
        erosion_large_aggregates,
    )
    deposition_gravel =
        natural_deposition(xs, median_diameter_gravel[v], input_gravel, erosion_gravel)

    return (
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_small_aggregates,
        deposition_large_aggregates,
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
    dt::Float64,
    v::Int,
)
    (;
        store_clay,
        store_silt,
        store_sand,
        store_small_aggregates,
        store_large_aggregates,
        store_gravel,
        deposition,
        erosion,
        sediment_rate,
        clay_rate,
        silt_rate,
        sand_rate,
        small_aggregates_rate,
        large_aggregates_rate,
        gravel_rate,
        leftover_clay,
        leftover_silt,
        leftover_sand,
        leftover_small_aggregates,
        leftover_large_aggregates,
        leftover_gravel,
    ) = variables

    (input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel) =
        input_particles
    (
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_small_aggregates,
        erosion_large_aggregates,
        erosion_gravel,
    ) = erosion_particles
    (
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_small_aggregates,
        deposition_large_aggregates,
        deposition_gravel,
    ) = deposition_particles

    # Update the sediment store
    store_clay[v] += deposition_clay * dt
    store_silt[v] += deposition_silt * dt
    store_sand[v] += deposition_sand * dt
    store_small_aggregates[v] += deposition_small_aggregates * dt
    store_large_aggregates[v] += deposition_large_aggregates * dt
    store_gravel[v] += deposition_gravel * dt

    # Compute total erosion
    erosion[v] = sum(erosion_particles)

    # Compute total deposition
    deposition[v] = sum(deposition_particles)

    # Output loads
    clay_rate[v] = fwaterout * (input_clay + erosion_clay - deposition_clay)
    silt_rate[v] = fwaterout * (input_silt + erosion_silt - deposition_silt)
    sand_rate[v] = fwaterout * (input_sand + erosion_sand - deposition_sand)
    small_aggregates_rate[v] =
        fwaterout * (input_sagg + erosion_small_aggregates - deposition_small_aggregates)
    large_aggregates_rate[v] =
        fwaterout * (input_lagg + erosion_large_aggregates - deposition_large_aggregates)
    gravel_rate[v] = fwaterout * (input_gravel + erosion_gravel - deposition_gravel)

    sediment_rate[v] =
        clay_rate[v] +
        silt_rate[v] +
        sand_rate[v] +
        small_aggregates_rate[v] +
        large_aggregates_rate[v] +
        gravel_rate[v]

    # Sediment left in the cell
    leftover_clay[v] = (input_clay + erosion_clay - deposition_clay - clay_rate[v]) * dt
    leftover_silt[v] = (input_silt + erosion_silt - deposition_silt - silt_rate[v]) * dt
    leftover_sand[v] = (input_sand + erosion_sand - deposition_sand - sand_rate[v]) * dt
    leftover_small_aggregates[v] =
        (
            input_sagg + erosion_small_aggregates - deposition_small_aggregates -
            small_aggregates_rate[v]
        ) * dt
    leftover_large_aggregates[v] =
        (
            input_lagg + erosion_large_aggregates - deposition_large_aggregates -
            large_aggregates_rate[v]
        ) * dt
    leftover_gravel[v] =
        (input_gravel + erosion_gravel - deposition_gravel - gravel_rate[v]) * dt

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
    (;
        store_clay,
        store_silt,
        store_sand,
        store_small_aggregates,
        store_large_aggregates,
        store_gravel,
    ) = sediment_transport_model.variables

    (; graph, order) = domain.network
    (; flow_width, flow_length, reservoir_coverage) = domain.parameters

    # Sediment transport - water balance in the river
    for v in order
        ### Sediment input in the cell (left from previous time step + from land + from upstream outflux) ###
        input_particles = compute_sediment_input(sediment_transport_model, graph, dt, v)
        input_sediment = sum(input_particles)

        ### River erosion ###
        # Erosion only if the load is below the transport capacity of the flow.
        sediment_need = if reservoir_coverage[v]
            # No erosion in reservoirs
            0.0
        else
            # Erosion only if the load is below the transport capacity of the flow.
            max(transport_capacity[v] - input_sediment, 0.0)
        end

        # Available sediment stored from previous deposition
        store_sediment =
            store_clay[v] +
            store_silt[v] +
            store_sand[v] +
            store_small_aggregates[v] +
            store_large_aggregates[v] +
            store_gravel[v]

        # Direct erosion from the river bed/bank
        erosion_particles = compute_direct_river_erosion(
            sediment_transport_model,
            sediment_need,
            store_sediment,
            dt,
            v,
        )

        # Erosion/degradation of the previously deposited sediment (from clay to gravel)
        store_erosion_particles =
            compute_store_erosion!(sediment_transport_model.variables, sediment_need, dt, v)

        # Update total erosion
        erosion_particles = erosion_particles .+ store_erosion_particles

        ### Deposition / settling ###

        # Different deposition if reservoir outlet or river
        deposition_particles = if reservoir_outlet[v]
            # Deposition in reservoir outlets
            compute_reservoir_deposition(
                sediment_transport_model,
                domain.parameters,
                input_particles,
                erosion_particles,
                v,
            )
        elseif reservoir_coverage[v]
            # No deposition in reservoir coverage, only at the outlets
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        else
            # Deposition in the river
            excess_sediment = max(input_sediment - transport_capacity[v], 0.0)

            if excess_sediment > 0.0
                # From transport capacity exceedance
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
                    v,
                )
            end
        end

        fwaterout =
            water_outflow_fraction(waterlevel[v], q[v], flow_width[v], flow_length[v], dt)

        update_variables!(
            sediment_transport_model.variables,
            input_particles,
            erosion_particles,
            deposition_particles,
            fwaterout,
            dt,
            v,
        )
    end
    return nothing
end

"Struct to store river sediment concentrations model variables"
@with_kw struct SedimentConcentrationsRiverVariables
    n_river_cells::Int
    # Total sediment concentration in the river [kg m⁻³]
    total::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # suspended sediment concentration in the river [kg m⁻³]
    suspended::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # bed load sediment concentration in the river [kg m⁻³]
    bed::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
end

"Struct to store river sediment concentrations model boundary conditions"
@with_kw struct SedimentConcentrationsRiverBC
    n_river_cells::Int
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n_river_cells) # [m]
    # Clay load [kg s⁻¹]
    clay::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Silt load [kg s⁻¹]
    silt::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Sand load [kg s⁻¹]
    sand::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Small aggregates load [kg s⁻¹]
    small_aggregates::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Large aggregates load [kg s⁻¹]
    large_aggregates::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Gravel load [kg s⁻¹]
    gravel::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
end

"Struct to store river sediment concentrations model parameters"
@with_kw struct SedimentConcentrationsRiverParameters
    # Clay mean diameter [m]
    median_diameter_clay::Vector{Float64}
    # Silt mean diameter [m]
    median_diameter_silt::Vector{Float64}
    # Sand mean diameter [m]
    median_diameter_sand::Vector{Float64}
    # Small aggregates mean diameter [m]
    median_diameter_small_aggregates::Vector{Float64}
    # Large aggregates mean diameter [m]
    median_diameter_large_aggregates::Vector{Float64}
    # Gravel mean diameter [m]
    median_diameter_gravel::Vector{Float64}
end

"Initialize river sediment concentrations model parameters"
function SedimentConcentrationsRiverParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    median_diameter_clay =
        ncread(dataset, config, "clay__mean_diameter", SoilLossModel; sel = indices)
    median_diameter_silt =
        ncread(dataset, config, "silt__mean_diameter", SoilLossModel; sel = indices)
    median_diameter_sand =
        ncread(dataset, config, "sand__mean_diameter", SoilLossModel; sel = indices)
    median_diameter_small_aggregates = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter",
        SoilLossModel;
        sel = indices,
    )
    median_diameter_large_aggregates = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter",
        SoilLossModel;
        sel = indices,
    )
    median_diameter_gravel =
        ncread(dataset, config, "gravel__mean_diameter", SoilLossModel; sel = indices)
    conc_parameters = SedimentConcentrationsRiverParameters(;
        median_diameter_clay,
        median_diameter_silt,
        median_diameter_sand,
        median_diameter_small_aggregates,
        median_diameter_large_aggregates,
        median_diameter_gravel,
    )

    return conc_parameters
end

"Struct to store river sediment concentrations model"
@with_kw struct SedimentConcentrationsRiverModel <: AbstractSedimentConcentrationsRiverModel
    n_river_cells::Int
    boundary_conditions::SedimentConcentrationsRiverBC =
        SedimentConcentrationsRiverBC(; n_river_cells)
    parameters::SedimentConcentrationsRiverParameters
    variables::SedimentConcentrationsRiverVariables =
        SedimentConcentrationsRiverVariables(; n_river_cells)
end

"Initialize river sediment concentrations model"
function SedimentConcentrationsRiverModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n_river_cells = length(indices)
    parameters = SedimentConcentrationsRiverParameters(dataset, config, indices)
    sediment_concentrations_model =
        SedimentConcentrationsRiverModel(; n_river_cells, parameters)
    return sediment_concentrations_model
end

"Update boundary conditions for river sediment concentrations model"
function update_bc_river_sediment_concentration_model!(
    sediment_concentrations_model::SedimentConcentrationsRiverModel,
    hydrological_forcing::HydrologicalForcing,
    sediment_transport_model::AbstractSedimentRiverTransportModel,
)
    (; q, waterlevel, clay, silt, sand, small_aggregates, large_aggregates, gravel) =
        sediment_concentrations_model.boundary_conditions
    # Hydrological forcing
    (; q_river, waterlevel_river) = hydrological_forcing
    @. q = q_river
    @. waterlevel = waterlevel_river
    # Sediment flux per particle
    @. clay = sediment_transport_model.variables.clay_rate
    @. silt = sediment_transport_model.variables.silt_rate
    @. sand = sediment_transport_model.variables.sand_rate
    @. small_aggregates = sediment_transport_model.variables.small_aggregates_rate
    @. large_aggregates = sediment_transport_model.variables.large_aggregates_rate
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
    (; q, waterlevel, clay, silt, sand, small_aggregates, large_aggregates, gravel) =
        sediment_transport_model.boundary_conditions
    (;
        median_diameter_clay,
        median_diameter_silt,
        median_diameter_sand,
        median_diameter_small_aggregates,
        median_diameter_large_aggregates,
        median_diameter_gravel,
    ) = sediment_transport_model.parameters
    (; total, suspended, bed) = sediment_transport_model.variables
    (; slope) = parameters

    for (river_idx, flow) in enumerate(q)
        if flow > 0
            # Differentiation of bed and suspended load using Rouse number for suspension
            # threshold diameter between bed load and mixed load using Rouse number
            common_term =
                0.41 * sqrt(
                    GRAVITATIONAL_ACCELERATION * waterlevel[river_idx] * slope[river_idx],
                ) / STOKES_FACTOR
            dbedf = 1e-3 * sqrt(2.5 * common_term)
            # # threshold diameter between suspended load and mixed load using Rouse number
            dsuspf = 1e-3 * sqrt(1.2 * common_term)

            # Rouse with diameter
            SSclay = suspended_solid(
                median_diameter_clay[river_idx],
                dsuspf,
                dbedf,
                clay[river_idx],
            )
            SSsilt = suspended_solid(
                median_diameter_silt[river_idx],
                dsuspf,
                dbedf,
                silt[river_idx],
            )
            SSsand = suspended_solid(
                median_diameter_sand[river_idx],
                dsuspf,
                dbedf,
                sand[river_idx],
            )
            SSsagg = suspended_solid(
                median_diameter_small_aggregates[river_idx],
                dsuspf,
                dbedf,
                small_aggregates[river_idx],
            )
            SSlagg = suspended_solid(
                median_diameter_large_aggregates[river_idx],
                dsuspf,
                dbedf,
                large_aggregates[river_idx],
            )
            SSgrav = suspended_solid(
                median_diameter_gravel[river_idx],
                dsuspf,
                dbedf,
                gravel[river_idx],
            )

            to_conc = inv(flow)
            total_ =
                clay[river_idx] +
                silt[river_idx] +
                small_aggregates[river_idx] +
                sand[river_idx] +
                large_aggregates[river_idx] +
                gravel[river_idx]
            total[river_idx] = total_ * to_conc

            SS = SSclay + SSsilt + SSsand + SSsagg + SSlagg + SSgrav
            suspended[river_idx] = SS * to_conc
            bed[river_idx] = (total_ - SS) * to_conc
        else
            suspended[river_idx] = 0.0
            bed[river_idx] = 0.0
            total[river_idx] = 0.0
        end
    end
end
