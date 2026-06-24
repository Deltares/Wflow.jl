abstract type AbstractIrrigationModel end

struct NoIrrigationPaddyModel <: AbstractIrrigationModel
    n::Int
end
struct NoIrrigationNonPaddyModel <: AbstractIrrigationModel
    n::Int
end
struct NoNonIrrigationDemandModel <: AbstractDemandModel
    n::Int
end
struct NoAllocationLandModel <: AbstractAllocationModel
    n::Int
end
struct NoAllocationRiverModel <: AbstractAllocationModel
    n::Int
end

"Struct to store non-irrigation water demand variables"
@with_kw struct NonIrrigationDemandVariables
    # return flow [m s⁻¹]
    returnflow::Vector{Float64}
    # return flow fraction [-]
    returnflow_fraction::Vector{Float64}
end

"Struct to store prescribed water demand variables"
@with_kw struct PrescribedDemand
    # gross water demand [m s⁻¹]
    demand_gross::Vector{Float64}
    # net water demand [m s⁻¹]
    demand_net::Vector{Float64}
end

"Non-irrigation water demand model"
@with_kw struct NonIrrigationDemandModel <: AbstractDemandModel
    demand::PrescribedDemand
    variables::NonIrrigationDemandVariables
end

# wrapper methods
get_demand_gross(demand_model::NonIrrigationDemandModel) = demand_model.demand.demand_gross
get_demand_gross(demand_model::NoNonIrrigationDemandModel) = Zeros(demand_model.n)

"Initialize non-irrigation water demand model for a water use `sector`"
function NonIrrigationDemandModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    sector::AbstractString,
)
    demand_gross = ncread(
        dataset,
        config,
        "$(sector)__gross_water_demand_volume_flux",
        LandHydrologySBM;
        sel = indices,
    )
    demand_net = ncread(
        dataset,
        config,
        "$(sector)__net_water_demand_volume_flux",
        LandHydrologySBM;
        sel = indices,
    )
    n = length(indices)
    returnflow_f = return_flow_fraction.(demand_gross, demand_net)

    demand = PrescribedDemand(; demand_gross, demand_net)
    vars = NonIrrigationDemandVariables(;
        returnflow_fraction = returnflow_f,
        returnflow = zeros(n),
    )
    non_irrigation_demand = NonIrrigationDemandModel(; demand, variables = vars)

    return non_irrigation_demand
end

"Struct to store non-paddy irrigation model variables"
@with_kw struct NonPaddyVariables
    n::Int
    demand_gross::Vector{Float64} = fill(MISSING_VALUE, n)   # irrigation gross demand [m s⁻¹]
end

"Struct to store non-paddy irrigation model parameters"
@with_kw struct NonPaddyParameters
    irrigation_efficiency::Vector{Float64}        # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{Float64}      # maximum irrigation rate [m s⁻¹]
    irrigation_areas::Vector{Bool}                # irrigation areas [-]
    irrigation_trigger::Vector{Bool}              # irrigation on or off [-]
end

"Non-paddy (other crops than flooded rice) irrigation model"
@with_kw struct NonPaddyModel <: AbstractIrrigationModel
    n::Int
    parameters::NonPaddyParameters
    variables::NonPaddyVariables = NonPaddyVariables(; n)
end

"Initialize non-paddy irrigation model"
function NonPaddyModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    efficiency = ncread(
        dataset,
        config,
        "irrigated_non_paddy__irrigation_efficiency",
        LandHydrologySBM;
        sel = indices,
    )
    areas = ncread(
        dataset,
        config,
        "irrigated_non_paddy_area__count",
        LandHydrologySBM;
        sel = indices,
    )
    irrigation_trigger = ncread(
        dataset,
        config,
        "irrigated_non_paddy__irrigation_trigger_flag",
        LandHydrologySBM;
        sel = indices,
    )
    max_irri_rate = ncread(
        dataset,
        config,
        "irrigated_non_paddy__max_irrigation_rate",
        LandHydrologySBM;
        sel = indices,
    )

    parameters = NonPaddyParameters(;
        maximum_irrigation_rate = max_irri_rate,
        irrigation_efficiency = efficiency,
        irrigation_areas = areas,
        irrigation_trigger,
    )
    n = length(indices)
    nonpaddy_model = NonPaddyModel(; n, parameters)

    return nonpaddy_model
end

# wrapper methods
get_demand_gross(nonpaddy_model::NonPaddyModel) = nonpaddy_model.variables.demand_gross
get_demand_gross(irrigation_model::NoIrrigationNonPaddyModel) = Zeros(irrigation_model.n)

"""
    update_demand_gross!(nonpaddy_model::NonPaddy, soil_model::SbmSoilModel, dt::Float64)

Update gross water demand `demand_gross` of the non-paddy irrigation model for a single
timestep.

The gross water demand is based on irrigation that is applied when the `irrigation_trigger`
is `true` (`on`) and when water depletion exceeds the readily available water in the root
zone of the SBM soil model. Irrigation brings the root zone back to field capacity, limited
by the infiltration capacity, taking into account limited irrigation efficiency and limited
by a maximum irrigation rate.
"""
function update_demand_gross!(
    nonpaddy_model::NonPaddyModel,
    soil_model::SbmSoilModel,
    dt::Float64,
)
    (; parameters, variables) = nonpaddy_model
    (; irrigation_areas, irrigation_trigger, maximum_irrigation_rate) = parameters
    (; demand_gross) = variables
    (; n_unsatlayers) = soil_model.variables

    for idx in eachindex(irrigation_areas)
        if irrigation_areas[idx] && irrigation_trigger[idx]
            irri_dem_gross_depth = 0.0
            for soil_layer_idx in 1:n_unsatlayers[idx]
                depletion, readily_available_water =
                    water_demand_root_zone(soil_model, idx, soil_layer_idx)

                # check if maximum irrigation rate has been applied at the previous time step.
                max_irri_rate_applied = demand_gross[idx] == maximum_irrigation_rate[idx]
                if depletion >= readily_available_water # start irrigation
                    irri_dem_gross_depth += depletion
                    # add depletion to irrigation gross demand when the maximum irrigation rate has been
                    # applied at the previous time step (to get volumetric water content at field capacity)
                elseif depletion > 0.0 && max_irri_rate_applied # continue irrigation
                    irri_dem_gross_depth += depletion
                end
            end
            irri_dem_gross = irri_dem_gross_depth / dt
            demand_gross[idx] =
                compute_demand_gross(nonpaddy_model, soil_model, irri_dem_gross, idx)
        else
            demand_gross[idx] = 0.0
        end
    end
    return nothing
end

update_demand_gross!(
    nonpaddy_model::NoIrrigationNonPaddyModel,
    soil_model::SbmSoilModel,
    dt::Float64,
) = nothing

"Compute water demand only for root zone through root fraction per layer"
function water_demand_root_zone(soil_model::SbmSoilModel, i::Int, k::Int)
    (;
        cumulative_layer_depth,
        air_entry_pressure,
        theta_s,
        theta_r,
        theta_fc,
        brooks_corey_exponent,
    ) = soil_model.parameters
    (; unsaturated_layer_thickness, unsaturated_layer_depth, h3) = soil_model.variables

    rooting_depth = get_rootingdepth(soil_model)

    rootfrac = min(
        1.0,
        (
            max(0.0, rooting_depth[i] - cumulative_layer_depth[i][k]) /
            unsaturated_layer_thickness[i][k]
        ),
    )
    # vwc_h3 can be precalculated.
    vwc_h3 = vwc_brooks_corey(
        h3[i],
        air_entry_pressure[i],
        theta_s[i],
        theta_r[i],
        brooks_corey_exponent[i][k],
    )
    depletion =
        (theta_fc[i] * unsaturated_layer_thickness[i][k]) -
        (unsaturated_layer_depth[i][k] + theta_r[i] * unsaturated_layer_thickness[i][k])
    depletion *= rootfrac
    readily_available_water = (theta_fc[i] - vwc_h3) * unsaturated_layer_thickness[i][k]
    readily_available_water *= rootfrac

    return depletion, readily_available_water
end

function compute_demand_gross(
    nonpaddy_model::NonPaddyModel,
    soil_model::SbmSoilModel,
    irri_dem_gross::Float64,
    i::Int,
)
    (; compacted_soil_area_fraction, infiltration_capacity_soil) = soil_model.parameters
    (; f_infiltration_reduction) = soil_model.variables
    (; irrigation_efficiency, maximum_irrigation_rate) = nonpaddy_model.parameters

    infiltration_capacity =
        f_infiltration_reduction[i] *
        (1.0 - compacted_soil_area_fraction[i]) *
        infiltration_capacity_soil[i]
    irri_dem_gross = min(irri_dem_gross, infiltration_capacity)
    irri_dem_gross /= irrigation_efficiency[i]
    # limit irrigation demand to the maximum irrigation rate
    irri_dem_gross = min(irri_dem_gross, maximum_irrigation_rate[i])

    return irri_dem_gross
end

"Struct to store paddy irrigation model variables"
@with_kw struct PaddyVariables
    n::Int
    demand_gross::Vector{Float64} = fill(MISSING_VALUE, n) # irrigation gross demand [m s⁻¹]
    h::Vector{Float64} = zeros(n)                          # actual water depth in rice field [m]
    evaporation::Vector{Float64} = zeros(n)                # evaporation rate [m s⁻¹]
end

"Struct to store paddy irrigation model parameters"
@with_kw struct PaddyParameters
    irrigation_efficiency::Vector{Float64}        # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{Float64}      # maximum irrigation rate [m s⁻¹]
    irrigation_areas::Vector{Bool}                # irrigation areas [-]
    irrigation_trigger::Vector{Bool}              # irrigation on or off [-]
    h_min::Vector{Float64}                        # minimum required water depth in the irrigated rice field [m]
    h_opt::Vector{Float64}                        # optimal water depth in the irrigated rice fields  [m]
    h_max::Vector{Float64}                        # water depth when rice field starts spilling water (overflow)  [m]
end

"PaddyModel (flooded rice) irrigation model"
@with_kw struct PaddyModel <: AbstractIrrigationModel
    n::Int
    parameters::PaddyParameters
    variables::PaddyVariables = PaddyVariables(; n)
end

"Initialize paddy irrigation model"
function PaddyModel(dataset::NCDataset, config::Config, indices::Vector{CartesianIndex{2}})
    h_min = ncread(
        dataset,
        config,
        "irrigated_paddy__min_depth",
        LandHydrologySBM;
        sel = indices,
    )
    h_opt = ncread(
        dataset,
        config,
        "irrigated_paddy__optimal_depth",
        LandHydrologySBM;
        sel = indices,
    )
    h_max = ncread(
        dataset,
        config,
        "irrigated_paddy__max_depth",
        LandHydrologySBM;
        sel = indices,
    )
    efficiency = ncread(
        dataset,
        config,
        "irrigated_paddy__irrigation_efficiency",
        LandHydrologySBM;
        sel = indices,
    )
    areas = ncread(
        dataset,
        config,
        "irrigated_paddy_area__count",
        LandHydrologySBM;
        sel = indices,
    )
    irrigation_trigger = ncread(
        dataset,
        config,
        "irrigated_paddy__irrigation_trigger_flag",
        LandHydrologySBM;
        sel = indices,
    )
    max_irri_rate = ncread(
        dataset,
        config,
        "irrigated_paddy__max_irrigation_rate",
        LandHydrologySBM;
        sel = indices,
    )
    n = length(indices)
    parameters = PaddyParameters(;
        irrigation_efficiency = efficiency,
        maximum_irrigation_rate = max_irri_rate,
        irrigation_trigger,
        h_min,
        h_max,
        h_opt,
        irrigation_areas = areas,
    )
    paddy = PaddyModel(; n, parameters)
    return paddy
end

# wrapper methods
get_water_depth(paddy_model::PaddyModel) = paddy_model.variables.h
get_water_depth(paddy_model::NoIrrigationPaddyModel) = Zeros(paddy_model.n)
get_demand_gross(paddy_model::PaddyModel) = paddy_model.variables.demand_gross
get_demand_gross(paddy_model::NoIrrigationPaddyModel) = Zeros(paddy_model.n)

"""
    evaporation!(paddy_model::Paddy, potential_evaporation, dt::Float64)

Update `evaporation` and the water depth `h` of the paddy irrigation model for a single
timestep.
"""
function evaporation!(paddy_model::PaddyModel, potential_evaporation, dt::Float64)
    for idx in eachindex(potential_evaporation)
        if paddy_model.parameters.irrigation_areas[idx]
            evaporation = min(paddy_model.variables.h[idx] / dt, potential_evaporation[idx])
            paddy_model.variables.h[idx] -= evaporation * dt
            paddy_model.variables.evaporation[idx] = evaporation
        end
    end
    return nothing
end
evaporation!(paddy_model::NoIrrigationPaddyModel, potential_evaporation, dt::Float64) =
    nothing

# wrapper methods
get_evaporation(paddy_model::NoIrrigationPaddyModel) = Zeros(paddy_model.n)
get_evaporation(paddy_model::PaddyModel) = paddy_model.variables.evaporation

"""
    update_runoff!(paddy_model::PaddyModel, runoff, dt::Float64)

Update `runoff` based on the water depth `h_max` (paddy field starts spilling), and update
the water depth `h` of the paddy irrigation model for a single timestep.
"""
function update_runoff!(paddy_model::PaddyModel, runoff, dt::Float64)
    (; parameters, variables) = paddy_model
    for (idx, b) in enumerate(parameters.irrigation_areas)
        !b && continue
        paddy_runoff = max(runoff[idx] - parameters.h_max[idx] / dt, 0.0)
        paddy_model.variables.h[idx] = (runoff[idx] - paddy_runoff) * dt
        runoff[idx] = paddy_runoff
    end
end

update_runoff!(paddy_model::NoIrrigationPaddyModel, runoff, dt::Float64) = nothing

"""
update_demand_gross!(paddy_model::Paddy, dt::Float64)

Update gross water demand `demand_gross` of the paddy irrigation model for a single
timestep.

The gross water demand is based on irrigation that is applied when the `irrigation_trigger`
is `true` (`on`) and when the paddy water depth `h` reaches below the minimum water depth
`h_min`. Irrigation is the amount required to reach the optimal paddy water depth `h_opt`,
taking into account limited irrigation efficiency and limited by a maximum irrigation rate.
"""
function update_demand_gross!(paddy_model::PaddyModel, dt::Float64)
    (; demand_gross) = paddy_model.variables
    (;
        irrigation_areas,
        irrigation_trigger,
        irrigation_efficiency,
        maximum_irrigation_rate,
    ) = paddy_model.parameters

    for idx in eachindex(irrigation_areas)
        if irrigation_areas[idx] && irrigation_trigger[idx]
            irr_depth_paddy = compute_irrigation_depth(paddy_model, idx)
            irri_dem_gross = irr_depth_paddy / irrigation_efficiency[idx]
            # limit irrigation demand to the maximum irrigation rate
            irri_dem_gross = min(irri_dem_gross / dt, maximum_irrigation_rate[idx])
            demand_gross[idx] = irri_dem_gross
        else
            demand_gross[idx] = 0.0
        end
    end
end

function compute_irrigation_depth(paddy_model::PaddyModel, i::Int)
    (; maximum_irrigation_rate, h_min, h_opt) = paddy_model.parameters
    (; demand_gross, h) = paddy_model.variables

    # check if maximum irrigation rate has been applied at the previous time step.
    max_irri_rate_applied = demand_gross[i] == maximum_irrigation_rate[i]
    # start irrigation
    irr_depth_paddy = if h[i] < h_min[i]
        h_opt[i] - h[i]
    elseif h[i] < h_opt[i] && max_irri_rate_applied # continue irrigation
        h_opt[i] - h[i]
    else
        0.0
    end

    return irr_depth_paddy
end

update_demand_gross!(paddy_model::NoIrrigationPaddyModel, dt::Float64) = nothing

"Struct to store water demand model variables"
@with_kw struct DemandVariables
    n::Int
    irrigation_demand_gross::Vector{Float64} = zeros(n)        # irrigation gross demand [m s⁻¹]
    non_irrigation_demand_gross::Vector{Float64} = zeros(n)     # non-irrigation gross demand [m s⁻¹]
    total_gross_demand::Vector{Float64} = zeros(n)       # total gross demand [m s¹]
    surfacewater_demand::Vector{Float64} = zeros(n)      # demand from surface water [m s⁻¹]
    groundwater_demand::Vector{Float64} = zeros(n)       # demand from groundwater [m s⁻¹]
end

"Water demand model"
@with_kw struct DemandModel{
    D <: AbstractDemandModel,
    I <: AbstractDemandModel,
    L <: AbstractDemandModel,
    P <: AbstractIrrigationModel,
    NP <: AbstractIrrigationModel,
} <: AbstractDemandModel
    domestic::D
    industry::I
    livestock::L
    paddy::P
    nonpaddy::NP
    variables::DemandVariables
end

@with_kw struct NoDemandModel <: AbstractDemandModel
    n::Int
    domestic::NoNonIrrigationDemandModel = NoNonIrrigationDemandModel(n)
    industry::NoNonIrrigationDemandModel = NoNonIrrigationDemandModel(n)
    livestock::NoNonIrrigationDemandModel = NoNonIrrigationDemandModel(n)
    paddy::NoIrrigationPaddyModel = NoIrrigationPaddyModel(n)
    nonpaddy::NoIrrigationNonPaddyModel = NoIrrigationNonPaddyModel(n)
end

"Initialize water demand model"
function DemandModel(dataset::NCDataset, config::Config, indices::Vector{CartesianIndex{2}})
    n = length(indices)
    demand(
        name;
        constr = NonIrrigationDemandModel,
        constr_triv = NoNonIrrigationDemandModel,
    ) =
        if getfield(config.model.water_demand, Symbol("$(name)__flag"))::Bool
            if constr == NonIrrigationDemandModel
                constr(dataset, config, indices, name)
            else
                constr(dataset, config, indices)
            end
        else
            constr_triv(n)
        end

    domestic = demand("domestic")
    industry = demand("industry")
    livestock = demand("livestock")
    paddy = demand("paddy"; constr = PaddyModel, constr_triv = NoIrrigationPaddyModel)
    nonpaddy =
        demand("nonpaddy"; constr = NonPaddyModel, constr_triv = NoIrrigationNonPaddyModel)

    variables = DemandVariables(; n)
    return DemandModel(; domestic, industry, livestock, paddy, nonpaddy, variables)
end

"Struct to store river allocation model variables"
@with_kw struct AllocationRiverVariables
    n_river::Int
    actual_surfacewater_abstraction::Vector{Float64} = zeros(n_river)       # actual surface water abstraction [m s⁻¹]
    actual_surfacewater_abstraction_volume::Vector{Float64} = zeros(n_river)   # actual surface water abstraction [m³ s⁻¹]
    available_surfacewater::Vector{Float64} = zeros(n_river)      # available surface water [m³]
    non_irrigation_returnflow::Vector{Float64} = zeros(n_river)          # return flow from non irrigation [m s⁻¹]
end

"River allocation model"
@with_kw struct AllocationRiverModel <: AbstractAllocationModel
    n_river::Int
    variables::AllocationRiverVariables = AllocationRiverVariables(; n_river)
end

get_nonirrigation_returnflow(allocation_model::AllocationRiverModel) =
    allocation_model.variables.non_irrigation_returnflow
get_nonirrigation_returnflow(allocation_model::NoAllocationRiverModel) =
    Zeros(allocation_model.n)

"Struct to store land allocation allocation model parameters"
@with_kw struct AllocationLandParameters
    fraction_surfacewater_used::Vector{Float64}     # fraction surface water used [-]
    areas::Vector{Int}          # allocation areas [-]
end

"Struct to store land allocation model variables"
@with_kw struct AllocationLandVariables
    n::Int
    # allocation from surface water [m s⁻¹]
    surfacewater_allocation::Vector{Float64} = zeros(n)
    # actual groundwater abstraction [m s⁻¹]
    actual_groundwater_abstraction::Vector{Float64} = zeros(n)
    # actual groundwater abstraction [m³ s⁻¹]
    actual_groundwater_abstraction_volume::Vector{Float64} = zeros(n)
    # available groundwater [m³]
    available_groundwater::Vector{Float64} = zeros(n)
    # allocation from groundwater [m s⁻¹]
    groundwater_allocation::Vector{Float64} = zeros(n)
    # allocated water for irrigation [m s⁻¹]
    irrigation_allocation::Vector{Float64} = zeros(n)
    # allocated water for non-irrigation [m s⁻¹]
    non_irrigation_allocation::Vector{Float64} = zeros(n)
    # total allocated water [m s⁻¹]
    total_alloc::Vector{Float64} = zeros(n)
    # return flow from non irrigation [m s⁻¹]
    non_irrigation_returnflow::Vector{Float64} = zeros(n)
end

"Land allocation model"
@with_kw struct AllocationLandModel <: AbstractAllocationModel
    n::Int
    parameters::AllocationLandParameters
    variables::AllocationLandVariables = AllocationLandVariables(; n)
end

"Initialize water allocation for the land domain"
function AllocationLandModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    fraction_surfacewater_used = ncread(
        dataset,
        config,
        "land_surface_water__withdrawal_fraction",
        LandHydrologySBM;
        sel = indices,
    )
    areas = ncread(
        dataset,
        config,
        "land_water_allocation_area__count",
        LandHydrologySBM;
        sel = indices,
    )

    # Check for allocation area ids that appear both inside and outside the active indices
    areas_full = ncread(
        dataset,
        config,
        "land_water_allocation_area__count",
        LandHydrologySBM;
        metadata = ParameterMetadata(; type = Int, allow_missing = true),
        logging = false,
    )
    ids_inside = Set(skipmissing(areas))
    mask_outside = trues(size(areas_full))
    mask_outside[indices] .= false
    ids_outside = Set(skipmissing(areas_full[mask_outside]))
    shared_ids = intersect(ids_inside, ids_outside)
    if !isempty(shared_ids)
        @warn "Allocation area ids `$(sort!(collect(shared_ids)))` are present both inside " *
              "and outside the active model domain. This may lead to incorrect water allocation."
    end

    n = length(indices)
    parameters = AllocationLandParameters(; areas, fraction_surfacewater_used)
    return AllocationLandModel(; n, parameters)
end

# wrapper methods
get_irrigation_allocated(allocation_model::AllocationLandModel) =
    allocation_model.variables.irrigation_allocation
get_irrigation_allocated(allocation_model::NoAllocationLandModel) =
    Zeros(allocation_model.n)
get_nonirrigation_returnflow(allocation_model::AllocationLandModel) =
    allocation_model.variables.non_irrigation_returnflow
get_nonirrigation_returnflow(allocation_model::NoAllocationLandModel) =
    Zeros(allocation_model.n)
get_groundwater_abstraction_flux(allocation_model::AllocationLandModel) =
    allocation_model.variables.actual_groundwater_abstraction
get_groundwater_abstraction_flux(allocation_model::NoAllocationLandModel) =
    Zeros(allocation_model.n)

"""
Return return flow fraction based on gross water demand `demand_gross` and net water demand
`demand_net`
"""
function return_flow_fraction(demand_gross::Float64, demand_net::Float64)
    fraction = bounded_divide(demand_net, demand_gross)
    returnflow_fraction = 1.0 - fraction
    return returnflow_fraction
end

"Update returnflow fraction for a non-irrigation water demand model"
function return_flow_fraction!(demand_model::NonIrrigationDemandModel)
    (; returnflow_fraction) = demand_model.variables
    (; demand_gross, demand_net) = demand_model.demand
    @. returnflow_fraction = return_flow_fraction(demand_gross, demand_net)
    return nothing
end

# return zero (gross water demand) if non-irrigation water demand sector is not defined
return_flow_fraction!(demand_model::NoNonIrrigationDemandModel) = nothing

"""
Update water allocation for river and land domains based on local surface water (river)
availability.
"""
function surface_water_allocation_local!(
    allocation_model::AllocationLandModel,
    demand_variables::DemandVariables,
    river_flow_model::AbstractRiverFlowModel,
    domain::DomainLand,
    dt::Float64,
)
    (; surfacewater_allocation) = allocation_model.variables
    (; surfacewater_demand) = demand_variables
    (;
        actual_surfacewater_abstraction_volume,
        actual_surfacewater_abstraction,
        available_surfacewater,
    ) = river_flow_model.allocation.variables
    (; external_inflow) = river_flow_model.boundary_conditions
    (; storage) = river_flow_model.variables
    (; area) = domain.parameters
    indices_river = domain.network.river_inds_excl_reservoir

    # maps from the land domain to the internal river domain (linear index), excluding reservoirs
    for (idx, river_idx) in enumerate(indices_river)
        if river_idx > 0
            # the available volume is limited by a fixed scaling factor of 0.8 to prevent
            # rivers completely drying out. check for abstraction through negative external
            # inflow first and adjust available volume.
            if external_inflow[river_idx] < 0.0
                available_volume = storage[river_idx] * 0.80
                max_river_abstraction =
                    min(-external_inflow[river_idx] * dt, available_volume)
                available_volume = max(available_volume - max_river_abstraction, 0.0)
            else
                available_volume = storage[river_idx] * 0.80
            end
            # satisfy surface water demand with available local river volume
            surfacewater_demand_vol = surfacewater_demand[idx] * area[idx]
            abstraction_vol = min(surfacewater_demand_vol * dt, available_volume)
            actual_surfacewater_abstraction_volume[river_idx] = abstraction_vol / dt
            # remaining available surface water and demand [m³]
            available_surfacewater[river_idx] = max(available_volume - abstraction_vol, 0.0)
            abstraction = abstraction_vol / (area[idx] * dt)
            surfacewater_demand[idx] = max(surfacewater_demand[idx] - abstraction, 0.0)
            # update actual abstraction from river and surface water allocation (land cell)
            actual_surfacewater_abstraction[river_idx] = abstraction
            surfacewater_allocation[idx] = abstraction
        end
    end
    return nothing
end

"""
Update water allocation for river and land domains based on surface water (river)
availability for allocation areas.
"""
function surface_water_allocation_area!(
    allocation_model::AllocationLandModel,
    demand_variables::DemandVariables,
    river_flow_model::AbstractRiverFlowModel,
    domain::Domain,
    dt::Float64,
)
    inds_river = domain.river.network.allocation_area_indices
    inds_land = domain.land.network.allocation_area_indices
    inds_reservoir = domain.river.network.reservoir_indices
    (; area) = domain.land.parameters

    (;
        available_surfacewater,
        actual_surfacewater_abstraction_volume,
        actual_surfacewater_abstraction,
    ) = river_flow_model.allocation.variables
    (; surfacewater_allocation) = allocation_model.variables
    (; surfacewater_demand) = demand_variables
    (; reservoir) = river_flow_model.boundary_conditions

    for alloc_area_idx in eachindex(inds_river)
        # surface water_demand (allocation area)
        sw_demand_vol = mapreduce(
            idx -> surfacewater_demand[idx] * area[idx],
            +,
            inds_land[alloc_area_idx],
        )

        sw_available = available_surface_water!(
            available_surfacewater,
            reservoir,
            inds_river[alloc_area_idx],
            inds_reservoir,
            dt,
        )

        # total actual surface water abstraction [m³] in an allocation area, minimum of
        # available surface water and demand in an allocation area.
        sw_abstraction = min(sw_available / dt, sw_demand_vol)

        # fraction of available surface water that can be abstracted at allocation area
        # level
        frac_abstract_sw = bounded_divide(sw_abstraction, sw_available / dt)
        # fraction of water demand that can be satisfied by available surface water at
        # allocation area level
        frac_allocate_sw = bounded_divide(sw_abstraction, sw_demand_vol)

        # water abstracted from surface water at each river cell (including reservoir
        # locations).
        for river_idx in inds_river[alloc_area_idx]
            actual_surfacewater_abstraction_volume[river_idx] +=
                frac_abstract_sw * available_surfacewater[river_idx] / dt
            actual_surfacewater_abstraction[river_idx] =
                actual_surfacewater_abstraction_volume[river_idx] /
                domain.river.parameters.cell_area[river_idx]
        end

        # water allocated to each land cell.
        for idx in inds_land[alloc_area_idx]
            surfacewater_allocation[idx] += frac_allocate_sw * surfacewater_demand[idx]
        end
    end
    return nothing
end

function available_surface_water!(
    available_surfacewater::Vector{Float64},
    reservoir_model,
    indices_river::Vector{Int},
    indices_reservoir::Vector{Int},
    dt::Float64,
)
    sw_available = 0.0
    for j in indices_river
        k = indices_reservoir[j]
        if k > 0
            # for reservoir locations use reservoir storage, check for abstraction
            # through external negative inflow first and adjust available volume
            external_inflow = reservoir_model.boundary_conditions.external_inflow[k]
            available_volume = reservoir_model.variables.storage[k] * 0.98
            if external_inflow < 0.0
                if available_volume > -external_inflow * dt
                    available_volume += external_inflow * dt
                else
                    available_volume = 0.0
                end
            end
            available_surfacewater[j] = available_volume
            sw_available += available_volume
        else
            # river volume
            sw_available += available_surfacewater[j]
        end
    end
    return sw_available
end

"Update water allocation for land domain based on local groundwater availability."
function groundwater_allocation_local!(
    allocation_model::AllocationLandModel,
    demand_variables::DemandVariables,
    groundwater_storage::Vector{Float64},
    parameters::LandParameters,
    dt,
)
    (;
        surfacewater_allocation,
        actual_groundwater_abstraction_volume,
        available_groundwater,
        actual_groundwater_abstraction,
        groundwater_allocation,
    ) = allocation_model.variables
    (; groundwater_demand, total_gross_demand) = demand_variables
    (; area, reservoir_coverage) = parameters

    for idx in eachindex(groundwater_demand)
        # groundwater demand based on allocation from surface water.

        groundwater_demand[idx] =
            max(total_gross_demand[idx] - surfacewater_allocation[idx], 0.0)
        # excluding reservoirs
        if !reservoir_coverage[idx]
            # satisfy groundwater demand with available local groundwater volume
            groundwater_demand_vol = groundwater_demand[idx] * area[idx]
            available_volume = groundwater_storage[idx] * 0.75 # limit available groundwater volume
            abstraction_vol = min(groundwater_demand_vol * dt, available_volume)
            gw_abst_vol_rate = abstraction_vol / dt
            actual_groundwater_abstraction_volume[idx] = gw_abst_vol_rate
            # remaining available groundwater and demand
            available_groundwater[idx] = max(available_volume - abstraction_vol, 0.0)
            abstraction = gw_abst_vol_rate / area[idx]
            groundwater_demand[idx] = max(groundwater_demand[idx] - abstraction, 0.0)
            # update actual abstraction from groundwater and groundwater allocation (land cell)
            actual_groundwater_abstraction[idx] = abstraction
            groundwater_allocation[idx] = abstraction
        end
    end
    return nothing
end

"""
Update water allocation for land domain based on groundwater availability for allocation
areas.

"""
function groundwater_allocation_area!(
    allocation_model::AllocationLandModel,
    demand_variables::DemandVariables,
    domain::Domain,
    dt::Float64,
)
    inds_river = domain.river.network.allocation_area_indices
    inds_land = domain.land.network.allocation_area_indices
    (;
        actual_groundwater_abstraction_volume,
        available_groundwater,
        actual_groundwater_abstraction,
        groundwater_allocation,
    ) = allocation_model.variables

    (; groundwater_demand) = demand_variables
    (; area) = domain.land.parameters

    # loop over allocation areas
    for alloc_area_idx in eachindex(inds_river)
        # groundwater demand and availability (allocation area)
        gw_demand_vol = 0.0
        gw_available = 0.0
        for idx in inds_land[alloc_area_idx]
            gw_demand_vol += groundwater_demand[idx] * area[idx]
            gw_available += available_groundwater[idx]
        end
        # total actual groundwater abstraction [m³] in an allocation area, minimum of
        # available  groundwater and demand in an allocation area.
        gw_abstraction = min(gw_available, gw_demand_vol * dt)

        # fraction of available groundwater that can be abstracted at allocation area level
        frac_abstract_gw = bounded_divide(gw_abstraction, gw_available)
        # fraction of water demand that can be satisfied by available groundwater at
        # allocation area level.
        frac_allocate_gw = bounded_divide(gw_abstraction, gw_demand_vol * dt)

        # water abstracted from groundwater and allocated.
        for idx in inds_land[alloc_area_idx]
            actual_groundwater_abstraction_volume[idx] +=
                frac_abstract_gw * available_groundwater[idx] / dt
            actual_groundwater_abstraction[idx] =
                actual_groundwater_abstraction_volume[idx] / area[idx]
            groundwater_allocation[idx] += frac_allocate_gw * groundwater_demand[idx]
        end
    end
    return nothing
end

"Return and update non-irrigation sector (domestic, livestock, industry) return flow"
function return_flow(
    demand_model::NonIrrigationDemandModel,
    non_irrigation_demand_gross::Vector{Float64},
    non_irrigation_allocation::Vector{Float64},
)
    for idx in eachindex(demand_model.variables.returnflow)
        frac = bounded_divide(
            demand_model.demand.demand_gross[idx],
            non_irrigation_demand_gross[idx],
        )
        allocate = frac * non_irrigation_allocation[idx]
        demand_model.variables.returnflow[idx] =
            demand_model.variables.returnflow_fraction[idx] * allocate
    end
    return demand_model.variables.returnflow
end

# return zero (return flow) if non-irrigation sector is not defined
return_flow(
    demand_model::NoNonIrrigationDemandModel,
    non_irrigation_demand_gross::Vector{Float64},
    non_irrigation_allocation::Vector{Float64},
) = 0.0

"""
    update_water_allocation_model!(
    allocation_model::AllocationLandModel,
    demand_model::DemandModel,
    routing::Routing,
    domain::Domain,
    dt::Float64,
)

Update water allocation for the land domain `AllocationLandModel` and water allocation for the
river domain (part of `routing`) based on the water `demand` model for a single timestep.
First, surface water abstraction is computed to satisfy local water demand (non-irrigation
and irrigation), and then updated (including reservoirs) to satisfy the remaining water
demand for allocation areas. Then groundwater abstraction is computed to satisfy the
remaining local water demand, and then updated to satisfy the remaining water demand for
allocation areas. Finally, non-irrigation return flows are updated.
"""
function update_water_allocation_model!(
    allocation_model::AllocationLandModel,
    demand_model::DemandModel,
    routing::Routing,
    domain::Domain,
    dt::Float64,
)
    river = routing.river_flow
    inds_river = domain.land.network.river_inds_excl_reservoir
    inds_reservoir = domain.reservoir.network.river_indices
    (;
        groundwater_allocation,
        surfacewater_allocation,
        actual_groundwater_abstraction,
        actual_groundwater_abstraction_volume,
        total_alloc,
        irrigation_allocation,
        non_irrigation_allocation,
        non_irrigation_returnflow,
    ) = allocation_model.variables

    (;
        surfacewater_demand,
        non_irrigation_demand_gross,
        irrigation_demand_gross,
        total_gross_demand,
    ) = demand_model.variables

    (; fraction_surfacewater_used) = allocation_model.parameters
    (; actual_surfacewater_abstraction, actual_surfacewater_abstraction_volume) =
        river.allocation.variables
    (; abstraction, reservoir) = river.boundary_conditions

    surfacewater_allocation .= 0.0
    actual_surfacewater_abstraction .= 0.0
    actual_surfacewater_abstraction_volume .= 0.0
    # total surface water demand for each land cell
    @. surfacewater_demand =
        fraction_surfacewater_used * (non_irrigation_demand_gross + irrigation_demand_gross)

    # local surface water demand and allocation (river, excluding reservoirs)
    surface_water_allocation_local!(
        allocation_model,
        demand_model.variables,
        river,
        domain.land,
        dt,
    )
    # surface water demand and allocation for areas
    surface_water_allocation_area!(
        allocation_model,
        demand_model.variables,
        river,
        domain,
        dt,
    )

    @. abstraction = actual_surfacewater_abstraction_volume

    # for reservoir locations set river abstraction at zero and abstract volume
    # from reservoir, including an update of waterlevel
    if !isnothing(reservoir)
        @. abstraction[inds_reservoir] = 0.0
        @. reservoir.variables.storage -=
            actual_surfacewater_abstraction_volume[inds_reservoir] * dt
        @. reservoir.variables.waterlevel = waterlevel(
            reservoir.parameters.storage_curve_type,
            reservoir.parameters.area,
            reservoir.variables.storage,
            reservoir.parameters.storage_waterlevel_curve,
        )
    end

    groundwater_allocation .= 0.0
    actual_groundwater_abstraction_volume .= 0.0
    actual_groundwater_abstraction .= 0.0
    # local groundwater demand and allocation
    groundwater_allocation_local!(
        allocation_model,
        demand_model.variables,
        routing.subsurface_flow.variables.storage,
        domain.land.parameters,
        dt,
    )
    # groundwater demand and allocation for areas
    groundwater_allocation_area!(allocation_model, demand_model.variables, domain, dt)

    # irrigation allocation
    for i in eachindex(total_alloc)
        total_alloc[i] = groundwater_allocation[i] + surfacewater_allocation[i]
        frac_irri = bounded_divide(irrigation_demand_gross[i], total_gross_demand[i])
        irrigation_allocation[i] = frac_irri * total_alloc[i]
        non_irrigation_allocation[i] = total_alloc[i] - irrigation_allocation[i]
    end

    # non-irrigation return flows
    returnflow_livestock = return_flow(
        demand_model.livestock,
        non_irrigation_demand_gross,
        non_irrigation_allocation,
    )
    returnflow_domestic = return_flow(
        demand_model.domestic,
        non_irrigation_demand_gross,
        non_irrigation_allocation,
    )
    returnflow_industry = return_flow(
        demand_model.industry,
        non_irrigation_demand_gross,
        non_irrigation_allocation,
    )

    @. non_irrigation_returnflow =
        returnflow_livestock + returnflow_domestic + returnflow_industry

    for i in eachindex(non_irrigation_returnflow)
        if inds_river[i] > 0.0
            k = inds_river[i]
            river.allocation.variables.non_irrigation_returnflow[k] =
                non_irrigation_returnflow[i]
            non_irrigation_returnflow[i] = 0.0
        end
    end
end

update_water_allocation_model!(
    allocation_model::NoAllocationLandModel,
    demand_model::NoDemandModel,
    routing::Routing,
    domain::Domain,
    dt::Float64,
) = nothing

"""
    update_demand_gross!(demand_model::DemandModel, dt::Float64)

Update total irrigation gross water demand `irrigation_demand_gross`, total non-irrigation gross
water demand `non_irrigation_demand_gross` and total gross water demand `total_gross_demand`.
"""
function update_demand_gross!(demand_model::DemandModel, dt::Float64)
    (; nonpaddy, paddy, domestic, industry, livestock) = demand_model
    (; irrigation_demand_gross, non_irrigation_demand_gross, total_gross_demand) =
        demand_model.variables
    # get gross water demands
    industry_dem = get_demand_gross(industry)
    domestic_dem = get_demand_gross(domestic)
    livestock_dem = get_demand_gross(livestock)
    nonpaddy_dem_gross = get_demand_gross(nonpaddy)
    paddy_dem_gross = get_demand_gross(paddy)
    # update gross water demands
    @. irrigation_demand_gross = nonpaddy_dem_gross + paddy_dem_gross
    @. non_irrigation_demand_gross = industry_dem + domestic_dem + livestock_dem
    @. total_gross_demand =
        nonpaddy_dem_gross + paddy_dem_gross + industry_dem + domestic_dem + livestock_dem

    return nothing
end

update_demand_gross!(demand_model::NoDemandModel, dt::Float64) = nothing

"""
    update_water_demand_model!(demand_model::DemandModel, soil_model::SbmSoilModel, dt::Float64)

Update the return flow fraction `returnflow_fraction` of `industry`, `domestic` and
`livestock`, gross water demand `demand_gross` of `paddy` and `nonpaddy` models, and the
total gross water demand, total irrigation gross water demand and total non-irrigation gross
water demand as part of the water `demand` model.
"""
function update_water_demand_model!(
    demand_model::DemandModel,
    soil_model::SbmSoilModel,
    dt::Float64,
)
    (; nonpaddy, paddy, domestic, industry, livestock) = demand_model

    return_flow_fraction!(industry)
    return_flow_fraction!(domestic)
    return_flow_fraction!(livestock)

    update_demand_gross!(nonpaddy, soil_model, dt)
    update_demand_gross!(paddy, dt)
    update_demand_gross!(demand_model, dt)

    return nothing
end
update_water_demand_model!(
    demand_model::NoDemandModel,
    soil_model::SbmSoilModel,
    dt::Float64,
) = nothing
