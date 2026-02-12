abstract type AbstractSedimentRiverTransportModel end
abstract type AbstractSedimentConcentrationsRiverModel end

"Struct to store river sediment transport model variables"
@with_kw struct SedimentRiverTransportVariables
    n::Int
    # Sediment flux [t dt⁻¹ => kg s⁻¹]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    clay_rate::Vector{Float64} = zeros(n)
    silt_rate::Vector{Float64} = zeros(n)
    sand_rate::Vector{Float64} = zeros(n)
    sagg_rate::Vector{Float64} = zeros(n)
    lagg_rate::Vector{Float64} = zeros(n)
    gravel_rate::Vector{Float64} = zeros(n)
    # Total Sediment deposition rate [t dt⁻¹ => kg s⁻¹]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total sediment erosion rate (from store + direct river bed/bank) [t dt⁻¹ => kg s⁻¹]
    erosion::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sediment / particle left in the cell [t => kg] - states
    leftover_clay::Vector{Float64} = zeros(n)
    leftover_silt::Vector{Float64} = zeros(n)
    leftover_sand::Vector{Float64} = zeros(n)
    leftover_sagg::Vector{Float64} = zeros(n)
    leftover_lagg::Vector{Float64} = zeros(n)
    leftover_gravel::Vector{Float64} = zeros(n)
    # Sediment / particle stored on the river bed after deposition [t => kg] -states
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
    # Transport capacity of the flow [t dt⁻¹ => kg s⁻¹]
    transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sediment input rate from land erosion [t dt⁻¹ => kg s⁻¹]
    erosion_land_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    erosion_land_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    erosion_land_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    erosion_land_sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    erosion_land_lagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sediment erosion rate from direct river erosion [t dt⁻¹ => kg s⁻¹]
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
    # Clay mean diameter [µm => m]
    dm_clay::Vector{Float64}
    # Silt mean diameter [µm => m]
    dm_silt::Vector{Float64}
    # Sand mean diameter [µm => m]
    dm_sand::Vector{Float64}
    # Small aggregates mean diameter [µm => m]
    dm_sagg::Vector{Float64}
    # Large aggregates mean diameter [µm => m]
    dm_lagg::Vector{Float64}
    # Gravel mean diameter [µm => m]
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
        "river_bottom_and_bank_clay__mass_fraction",
        SoilLoss;
        sel=indices,
    )
    silt_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_silt__mass_fraction",
        SoilLoss;
        sel=indices,
    )
    sand_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_sand__mass_fraction",
        SoilLoss;
        sel=indices,
    )
    gravel_fraction = ncread(
        dataset,
        config,
        "river_bottom_and_bank_gravel__mass_fraction",
        SoilLoss;
        sel=indices,
    )
    # Check that river fractions sum to 1
    river_fractions = clay_fraction + silt_fraction + sand_fraction + gravel_fraction
    if any(abs.(river_fractions .- 1.0) .> 1e-3)
        error("Particle fractions in the river bed must sum to 1")
    end
    dm_clay = ncread(dataset, config, "clay__mean_diameter", SoilLoss; sel=indices)
    dm_silt = ncread(dataset, config, "silt__mean_diameter", SoilLoss; sel=indices)
    dm_sand = ncread(dataset, config, "sand__mean_diameter", SoilLoss; sel=indices)
    dm_sagg = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter",
        SoilLoss;
        sel=indices,
    )
    dm_lagg = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter",
        SoilLoss;
        sel=indices,
    )
    dm_gravel = ncread(dataset, config, "gravel__mean_diameter", SoilLoss; sel=indices)

    # Reservoirs
    if config.model.reservoir__flag
        reservoir_outlet =
            ncread(dataset, config, "reservoir_location__count", SoilLoss; sel=indices)
        reservoir_area =
            ncread(dataset, config, "reservoir_surface__area", SoilLoss; sel=indices)
        reservoir_trapping_efficiency = ncread(
            dataset,
            config,
            "reservoir_water_sediment__bedload_trapping_efficiency",
            SoilLoss;
            sel=indices,
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
        reservoir_outlet=reservoir_outlet .> 0,
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
    potential_erosion_model::RiverErosionJulianTorresModel,
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
    # [m³ s⁻¹] = [m³ s⁻¹]
    @. q = q_river
    # [m] = [m]
    @. waterlevel = waterlevel_river
    # Transport capacity
    # [kg s⁻¹] = [kg s⁻¹]
    @. transport_capacity = transport_capacity_model.variables.sediment_transport_capacity
    # Input from soil erosion
    (; clay_rate, silt_rate, sand_rate, sagg_rate, lagg_rate) = to_river_model.variables
    # [kg s⁻¹] = [kg s⁻¹]
    map!(i -> clay_rate[i], erosion_land_clay, indices_riv)
    map!(i -> silt_rate[i], erosion_land_silt, indices_riv)
    map!(i -> sand_rate[i], erosion_land_sand, indices_riv)
    map!(i -> sagg_rate[i], erosion_land_sagg, indices_riv)
    map!(i -> lagg_rate[i], erosion_land_lagg, indices_riv)
    # Maximum direct river bed/bank erosion
    # [kg s⁻¹] = [kg s⁻¹]
    @. potential_erosion_river_bed = potential_erosion_model.variables.bed
    @. potential_erosion_river_bank = potential_erosion_model.variables.bank
end

"""
Calculate sediment input from leftover sediment, land erosion, and upstream contributions
"""
function compute_sediment_input(
    model::SedimentRiverTransportModel,
    graph::DiGraph,
    dt::Float64,
    v::Int,
)
    (; boundary_conditions, variables) = model
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
    # [kg s⁻¹] = [kg] / [s] + [kg s⁻¹]
    input_clay = leftover_clay[v] / dt + erosion_land_clay[v]
    input_silt = leftover_silt[v] / dt + erosion_land_silt[v]
    input_sand = leftover_sand[v] / dt + erosion_land_sand[v]
    input_sagg = leftover_sagg[v] / dt + erosion_land_sagg[v]
    input_lagg = leftover_lagg[v] / dt + erosion_land_lagg[v]
    input_gravel = leftover_gravel[v] / dt

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

    return (input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel)
end

"""
Calculate direct river bed/bank erosion based on sediment need
"""
function compute_direct_river_erosion(
    model::SedimentRiverTransportModel,
    sediment_need::Float64,
    store_sediment::Float64,
    dt::Float64,
    v::Int,
)
    (; potential_erosion_river_bed, potential_erosion_river_bank) =
        model.boundary_conditions
    (; clay_fraction, silt_fraction, sand_fraction, gravel_fraction) = model.parameters

    if sediment_need <= store_sediment / dt
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    # Effective sediment needed from river bed and bank erosion
    # [kg s⁻¹] = [kg s⁻¹] - [kg] / [s]
    effsediment_need = sediment_need - store_sediment / dt

    # Relative potential erosion rates of the bed and the bank
    # [-]
    RTEbank = if (potential_erosion_river_bank[v] + potential_erosion_river_bed[v] > 0.0)
        # [kg s⁻¹] / ([kg s⁻¹] + [kg s⁻¹])
        potential_erosion_river_bank[v] /
        (potential_erosion_river_bank[v] + potential_erosion_river_bed[v])
    else
        0.0
    end
    # [-]
    RTEbed = 1.0 - RTEbank

    # Actual bed and bank erosion
    # [kg s⁻¹] = min([-] * [kg s⁻¹], [kg s⁻¹])
    erosion_bank = min(RTEbank * effsediment_need, potential_erosion_river_bank[v])
    erosion_bed = min(RTEbed * effsediment_need, potential_erosion_river_bed[v])
    # [kg s⁻¹] = [kg s⁻¹] + [kg s⁻¹]
    erosion_river = erosion_bank + erosion_bed
    # Per particle
    # [kg s⁻¹] = [kg s⁻¹] * [-]
    erosion_clay = erosion_river * clay_fraction[v]
    erosion_silt = erosion_river * silt_fraction[v]
    erosion_sand = erosion_river * sand_fraction[v]
    erosion_gravel = erosion_river * gravel_fraction[v]
    # No small and large aggregates in the river bed/bank
    # [kg s⁻¹]
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
    dt::Float64,
    v::Int,
)
    # [kg]
    (; store_clay, store_silt, store_sand, store_sagg, store_lagg, store_gravel) = variables

    # [kg s⁻¹], [kg s⁻¹]
    erosion_clay, sediment_need = river_erosion_store!(store_clay, sediment_need, dt, v)
    erosion_silt, sediment_need = river_erosion_store!(store_silt, sediment_need, dt, v)
    erosion_sagg, sediment_need = river_erosion_store!(store_sagg, sediment_need, dt, v)
    erosion_sand, sediment_need = river_erosion_store!(store_sand, sediment_need, dt, v)
    erosion_lagg, sediment_need = river_erosion_store!(store_lagg, sediment_need, dt, v)
    erosion_gravel, _ = river_erosion_store!(store_gravel, sediment_need, dt, v)

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
    model::SedimentRiverTransportModel,
    domain_parameters::RiverParameters,
    input_particles::NTuple{6,Float64},
    erosion_particles::NTuple{6,Float64},
    v::Int,
)
    (; boundary_conditions, parameters) = model
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
        q[v],
        waterlevel[v],
        reservoir_area[v],
        reservoir_trapping_efficiency[v],
        dm[v],
        slope[v],
    )

    # [kg s⁻¹]
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
    # [kg s⁻¹] = [kg s⁻¹] + [kg s⁻¹]
    total = input + erosion

    if excess_sediment > total
        # [kg s⁻¹]
        deposition = total
        # [kg s⁻¹] -= [kg s⁻¹]
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
    # [kg s⁻¹]
    input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel =
        input_particles
    # [kg s⁻¹]
    erosion_clay, erosion_silt, erosion_sand, erosion_sagg, erosion_lagg, erosion_gravel =
        erosion_particles

    # From largest to smallest particles
    # [kg s⁻¹], [kg s⁻¹]
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
    # [-]
    x_material = min(1.0, 1.0 - exp(-xs * fall_velocity(dm)))
    # [-] * ([kg s⁻¹] + [kg s⁻¹])
    return x_material * (input + erosion)
end

"""
Calculate natural river deposition using Einstein's formula (Stokes settling)
"""
function compute_natural_deposition(
    model::SedimentRiverTransportModel,
    domain_parameters::RiverParameters,
    input_particles::NTuple{6,Float64},
    erosion_particles::NTuple{6,Float64},
    v::Int,
)
    (; boundary_conditions, parameters) = model
    (; q) = boundary_conditions
    (; dm_clay, dm_silt, dm_sand, dm_sagg, dm_lagg, dm_gravel) = parameters
    (; flow_length, flow_width) = domain_parameters

    input_clay, input_silt, input_sand, input_sagg, input_lagg, input_gravel =
        input_particles
    erosion_clay, erosion_silt, erosion_sand, erosion_sagg, erosion_lagg, erosion_gravel =
        erosion_particles

    # Particle fall velocity [m/s] from Stokes
    xs = ifelse(q[v] > 0.0, 1.055 * flow_length[v] / (q[v] / flow_width[v]), 0.0)

    # [kg s⁻¹]
    deposition_clay = natural_deposition(xs, dm_clay[v], input_clay, erosion_clay)
    deposition_silt = natural_deposition(xs, dm_silt[v], input_silt, erosion_silt)
    deposition_sand = natural_deposition(xs, dm_sand[v], input_sand, erosion_sand)
    deposition_sagg = natural_deposition(xs, dm_sagg[v], input_sagg, erosion_sagg)
    deposition_lagg = natural_deposition(xs, dm_lagg[v], input_lagg, erosion_lagg)
    deposition_gravel = natural_deposition(xs, dm_gravel[v], input_gravel, erosion_gravel)

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
        # [-] = min([m³ s⁻¹] * [s] / ([m] * [m] * [m]), [-])
        min(q * dt / (waterlevel * flow_width * flow_length), 1.0)
    else
        1.0
    end
end

function update_variables!(
    variables::SedimentRiverTransportVariables,
    input_particles::NTuple{6,Float64},
    erosion_particles::NTuple{6,Float64},
    deposition_particles::NTuple{6,Float64},
    fwaterout::Float64,
    dt::Float64,
    v::Int,
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
    # [kg] += [kg s⁻¹] * [s]
    store_clay[v] += deposition_clay * dt
    store_silt[v] += deposition_silt * dt
    store_sand[v] += deposition_sand * dt
    store_sagg[v] += deposition_sagg * dt
    store_lagg[v] += deposition_lagg * dt
    store_gravel[v] += deposition_gravel * dt

    # Compute total erosion
    # [kg s⁻¹] = ∑ [kg s⁻¹]
    erosion[v] = sum(erosion_particles)

    # Compute total deposition
    # [kg s⁻¹] = ∑ [kg s⁻¹]
    deposition[v] = sum(deposition_particles)

    # Output loads
    # [kg s⁻¹] = [-] * ([kg s⁻¹] + [kg s⁻¹] - [kg s⁻¹])
    clay_rate[v] = fwaterout * (input_clay + erosion_clay - deposition_clay)
    silt_rate[v] = fwaterout * (input_silt + erosion_silt - deposition_silt)
    sand_rate[v] = fwaterout * (input_sand + erosion_sand - deposition_sand)
    sagg_rate[v] = fwaterout * (input_sagg + erosion_sagg - deposition_sagg)
    lagg_rate[v] = fwaterout * (input_lagg + erosion_lagg - deposition_lagg)
    gravel_rate[v] = fwaterout * (input_gravel + erosion_gravel - deposition_gravel)

    # [kg s⁻¹] = ∑ [kg s⁻¹]
    sediment_rate[v] =
        clay_rate[v] +
        silt_rate[v] +
        sand_rate[v] +
        sagg_rate[v] +
        lagg_rate[v] +
        gravel_rate[v]

    # Sediment left in the cell
    # [kg] = ([kg s⁻¹] + [kg s⁻¹] - [kg s⁻¹] - [kg s⁻¹]) * dt
    leftover_clay[v] = (input_clay + erosion_clay - deposition_clay - clay_rate[v]) * dt
    leftover_silt[v] = (input_silt + erosion_silt - deposition_silt - silt_rate[v]) * dt
    leftover_sand[v] = (input_sand + erosion_sand - deposition_sand - sand_rate[v]) * dt
    leftover_sagg[v] = (input_sagg + erosion_sagg - deposition_sagg - sagg_rate[v]) * dt
    leftover_lagg[v] = (input_lagg + erosion_lagg - deposition_lagg - lagg_rate[v]) * dt
    leftover_gravel[v] =
        (input_gravel + erosion_gravel - deposition_gravel - gravel_rate[v]) * dt

    return nothing
end

"update river sediment transport model for a single timestep"
function update!(model::SedimentRiverTransportModel, domain::DomainRiver, dt::Float64)
    (; waterlevel, q, transport_capacity) = model.boundary_conditions
    (; reservoir_outlet) = model.parameters
    (; store_clay, store_silt, store_sand, store_sagg, store_lagg, store_gravel) =
        model.variables

    (; graph, order) = domain.network
    (; flow_width, flow_length, reservoir_coverage) = domain.parameters

    # Sediment transport - water balance in the river
    for v in order
        ### Sediment input in the cell (left from prevoous time step + from land + from upstream outflux) ###
        # [kg s⁻¹]
        input_particles = compute_sediment_input(model, graph, dt, v)
        # [kg s⁻¹] = ∑ [kg s⁻¹]
        input_sediment = sum(input_particles)

        ### River erosion ###
        # Erosion only if the load is below the transport capacity of the flow.
        # [kg s⁻¹]
        sediment_need = if reservoir_coverage[v]
            # No erosion in reservoirs
            0.0
        else
            # Erosion only if the load is below the transport capacity of the flow.
            # [kg s⁻¹] = max([kg s⁻¹] - [kg s⁻¹], [kg s⁻¹])
            max(transport_capacity[v] - input_sediment, 0.0)
        end

        # Available sediment stored from previous deposition
        # [kg] = ∑ [kg]
        store_sediment =
            store_clay[v] +
            store_silt[v] +
            store_sand[v] +
            store_sagg[v] +
            store_lagg[v] +
            store_gravel[v]

        # Direct erosion from the river bed/bank
        # [kg s⁻¹]
        erosion_particles =
            compute_direct_river_erosion(model, sediment_need, store_sediment, dt, v)

        # Erosion/degradation of the previously deposited sediment (from clay to gravel)
        # [kg s⁻¹]
        store_erosion_particles =
            compute_store_erosion!(model.variables, sediment_need, dt, v)

        # Update total erosion
        # [kg s⁻¹] = [kg s⁻¹] + [kg s⁻¹]
        erosion_particles = erosion_particles .+ store_erosion_particles

        ### Deposition / settling ###

        # Different deposition if reservoir outlet or river
        # [kg s⁻¹]
        deposition_particles = if reservoir_outlet[v]
            # Deposition in reservoir outlets
            compute_reservoir_deposition(
                model,
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
            # [kg s⁻¹] = max([kg s⁻¹] - [kg s⁻¹], [kg s⁻¹])
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
                    model,
                    domain.parameters,
                    input_particles,
                    erosion_particles,
                    v,
                )
            end
        end

        # [-]
        fwaterout =
            water_outflow_fraction(waterlevel[v], q[v], flow_width[v], flow_length[v], dt)

        update_variables!(
            model.variables,
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
    n::Int
    # Total sediment concentration in the river [g m⁻³ => kg m⁻³]
    total::Vector{Float64} = fill(MISSING_VALUE, n)
    # suspended sediment concentration in the river [g m⁻³ => kg m⁻³]
    suspended::Vector{Float64} = fill(MISSING_VALUE, n)
    # bed load sediment concentration in the river [g m⁻³ => kg m⁻³]
    bed::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store river sediment concentrations model boundary conditions"
@with_kw struct SedimentConcentrationsRiverBC
    n::Int
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n)
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n) # [m]
    # Clay load [g m⁻³ => kg m⁻³]
    clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt load [g m⁻³ => kg m⁻³]
    silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand load [g m⁻³ => kg m⁻³]
    sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates load [g m⁻³ => kg m⁻³]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates load [g m⁻³ => kg m⁻³]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Gravel load [g m⁻³ => kg m⁻³]
    gravel::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store river sediment concentrations model parameters"
@with_kw struct SedimentConcentrationsRiverParameters
    # Clay mean diameter [µm => m]
    dm_clay::Vector{Float64}
    # Silt mean diameter [µm => m]
    dm_silt::Vector{Float64}
    # Sand mean diameter [µm => m]
    dm_sand::Vector{Float64}
    # Small aggregates mean diameter [µm => m]
    dm_sagg::Vector{Float64}
    # Large aggregates mean diameter [µm => m]
    dm_lagg::Vector{Float64}
    # Gravel mean diameter [µm => m]
    dm_gravel::Vector{Float64}
end

"Initialize river sediment concentrations model parameters"
function SedimentConcentrationsRiverParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    dm_clay = ncread(dataset, config, "clay__mean_diameter", SoilLoss; sel=indices)
    dm_silt = ncread(dataset, config, "silt__mean_diameter", SoilLoss; sel=indices)
    dm_sand = ncread(dataset, config, "sand__mean_diameter", SoilLoss; sel=indices)
    dm_sagg = ncread(
        dataset,
        config,
        "sediment_small_aggregates__mean_diameter",
        SoilLoss;
        sel=indices,
    )
    dm_lagg = ncread(
        dataset,
        config,
        "sediment_large_aggregates__mean_diameter",
        SoilLoss;
        sel=indices,
    )
    dm_gravel = ncread(dataset, config, "gravel__mean_diameter", SoilLoss; sel=indices)
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
    sediment_flux_model::SedimentRiverTransportModel,
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
            # [m]
            dbedf = 1e-3 * sqrt(2.5 * common_term)
            # # threshold diameter between suspended load and mixed load using Rouse number
            # [m]
            dsuspf = 1e-3 * sqrt(1.2 * common_term)

            # Rouse with diameter
            # [kg m⁻³]
            SSclay = suspended_solid(dm_clay[i], dsuspf, dbedf, clay[i])
            SSsilt = suspended_solid(dm_silt[i], dsuspf, dbedf, silt[i])
            SSsand = suspended_solid(dm_sand[i], dsuspf, dbedf, sand[i])
            SSsagg = suspended_solid(dm_sagg[i], dsuspf, dbedf, sagg[i])
            SSlagg = suspended_solid(dm_lagg[i], dsuspf, dbedf, lagg[i])
            SSgrav = suspended_solid(dm_gravel[i], dsuspf, dbedf, gravel[i])

            to_conc = 1e6 / (flow * dt)
            # [kg m⁻³] = ∑ [kg m⁻³]
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
