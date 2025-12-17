abstract type AbstractIrrigationModel end

struct NoIrrigationPaddy <: AbstractIrrigationModel
    n::Int
end
struct NoIrrigationNonPaddy <: AbstractIrrigationModel
    n::Int
end
struct NoNonIrrigationDemand <: AbstractDemandModel
    n::Int
end
struct NoAllocationLand <: AbstractAllocationModel
    n::Int
end
struct NoAllocationRiver <: AbstractAllocationModel
    n::Int
end

"Struct to store non-irrigation water demand variables"
@with_kw struct NonIrrigationDemandVariables
    returnflow::Vector{Float64}               # return flow [mm Δt⁻¹]
    returnflow_fraction::Vector{Float64}      # return flow fraction [-]
end

"Struct to store prescribed water demand variables"
@with_kw struct PrescibedDemand
    demand_gross::Vector{Float64}     # gross water demand [mm Δt⁻¹]
    demand_net::Vector{Float64}       # net water demand [mm Δt⁻¹]
end

"Non-irrigation water demand model"
@with_kw struct NonIrrigationDemand <: AbstractDemandModel
    demand::PrescibedDemand
    variables::NonIrrigationDemandVariables
end

# wrapper methods
get_demand_gross(model::NonIrrigationDemand) = model.demand.demand_gross
get_demand_gross(model::NoNonIrrigationDemand) = Zeros(model.n)

"Initialize non-irrigation water demand model for a water use `sector`"
function NonIrrigationDemand(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
    sector::AbstractString,
)
    demand_gross =
        ncread(
            dataset,
            config,
            "$(sector)__gross_water_demand_volume_flux";
            sel = indices,
            defaults = 0.0,
            type = Float64,
        ) .* (dt / BASETIMESTEP)
    demand_net =
        ncread(
            dataset,
            config,
            "$(sector)__net_water_demand_volume_flux";
            sel = indices,
            defaults = 0.0,
            type = Float64,
        ) .* (dt / BASETIMESTEP)
    n = length(indices)
    returnflow_f = return_flow_fraction.(demand_gross, demand_net)

    demand = PrescibedDemand(; demand_gross, demand_net)
    vars = NonIrrigationDemandVariables(;
        returnflow_fraction = returnflow_f,
        returnflow = fill(Float64(0), n),
    )
    non_irrigation_demand = NonIrrigationDemand(; demand, variables = vars)

    return non_irrigation_demand
end

"Struct to store non-paddy irrigation model variables"
@with_kw struct NonPaddyVariables
    n::Int
    demand_gross::Vector{Float64} = fill(MISSING_VALUE, n)     # irrigation gross demand [mm Δt⁻¹]
end

"Struct to store non-paddy irrigation model parameters"
@with_kw struct NonPaddyParameters
    irrigation_efficiency::Vector{Float64}        # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{Float64}      # maximum irrigation rate [mm Δt⁻¹]
    irrigation_areas::Vector{Bool}          # irrigation areas [-]
    irrigation_trigger::Vector{Bool}        # irrigation on or off [-]
end

"Non-paddy (other crops than flooded rice) irrigation model"
@with_kw struct NonPaddy <: AbstractIrrigationModel
    parameters::NonPaddyParameters
    variables::NonPaddyVariables
end

"Initialize non-paddy irrigation model"
function NonPaddy(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    efficiency = ncread(
        dataset,
        config,
        "irrigated_non_paddy__irrigation_efficiency";
        sel = indices,
        defaults = 1.0,
        type = Float64,
    )
    areas = ncread(
        dataset,
        config,
        "irrigated_non_paddy_area__count";
        optional = false,
        sel = indices,
        type = Int,
    )
    irrigation_trigger = ncread(
        dataset,
        config,
        "irrigated_non_paddy__irrigation_trigger_flag";
        optional = false,
        sel = indices,
        type = Bool,
    )
    max_irri_rate =
        ncread(
            dataset,
            config,
            "irrigated_non_paddy__max_irrigation_rate";
            sel = indices,
            defaults = 25.0,
            type = Float64,
        ) .* (dt / BASETIMESTEP)

    parameters = NonPaddyParameters(;
        maximum_irrigation_rate = max_irri_rate,
        irrigation_efficiency = efficiency,
        irrigation_areas = areas,
        irrigation_trigger,
    )
    n = length(indices)
    variables = NonPaddyVariables(; n)

    nonpaddy = NonPaddy(; variables, parameters)

    return nonpaddy
end

# wrapper methods
get_demand_gross(model::NonPaddy) = model.variables.demand_gross
get_demand_gross(model::NoIrrigationNonPaddy) = Zeros(model.n)

"""
    update_demand_gross!(model::NonPaddy, soil::SbmSoilModel)

Update gross water demand `demand_gross` of the non-paddy irrigation model for a single
timestep.

The gross water demand is based on irrigation that is applied when the `irrigation_trigger`
is `true` (`on`) and when water depletion exceeds the readily available water in the root
zone of the SBM soil model. Irrigation brings the root zone back to field capacity, limited
by the infiltration capacity, taking into account limited irrigation efficiency and limited
by a maximum irrigation rate.
"""
function update_demand_gross!(model::NonPaddy, soil::SbmSoilModel)
    (; hb, theta_s, theta_r, c, sumlayers, act_thickl, pathfrac, infiltcapsoil) =
        soil.parameters
    (;
        h3,
        ustorelayerthickness,
        n_unsatlayers,
        ustorelayerdepth,
        f_infiltration_reduction,
    ) = soil.variables
    (;
        irrigation_areas,
        irrigation_trigger,
        maximum_irrigation_rate,
        irrigation_efficiency,
    ) = model.parameters
    rootingdepth = get_rootingdepth(soil)

    for i in eachindex(irrigation_areas)
        if irrigation_areas[i] && irrigation_trigger[i]
            irri_dem_gross = 0.0
            for k in 1:n_unsatlayers[i]
                # compute water demand only for root zone through root fraction per layer
                rootfrac = min(
                    1.0,
                    (
                        max(0.0, rootingdepth[i] - sumlayers[i][k]) /
                        ustorelayerthickness[i][k]
                    ),
                )
                # vwc_f and vwc_h3 can be precalculated.
                vwc_fc = vwc_brooks_corey(-100.0, hb[i], theta_s[i], theta_r[i], c[i][k])
                vwc_h3 = vwc_brooks_corey(h3[i], hb[i], theta_s[i], theta_r[i], c[i][k])
                depletion =
                    (vwc_fc * ustorelayerthickness[i][k]) -
                    (ustorelayerdepth[i][k] + theta_r[i] * ustorelayerthickness[i][k])
                depletion *= rootfrac
                raw = (vwc_fc - vwc_h3) * ustorelayerthickness[i][k] # readily available water
                raw *= rootfrac

                # check if maximum irrigation rate has been applied at the previous time step.
                max_irri_rate_applied =
                    model.variables.demand_gross[i] == maximum_irrigation_rate[i]
                if depletion >= raw # start irrigation
                    irri_dem_gross += depletion
                    # add depletion to irrigation gross demand when the maximum irrigation rate has been
                    # applied at the previous time step (to get volumetric water content at field capacity)
                elseif depletion > 0.0 && max_irri_rate_applied # continue irrigation
                    irri_dem_gross += depletion
                end
            end
            # limit irrigation demand to infiltration capacity
            infiltration_capacity =
                f_infiltration_reduction[i] * (1.0 - pathfrac[i]) * infiltcapsoil[i]
            irri_dem_gross = min(irri_dem_gross, infiltration_capacity)
            irri_dem_gross /= irrigation_efficiency[i]
            # limit irrigation demand to the maximum irrigation rate
            irri_dem_gross = min(irri_dem_gross, maximum_irrigation_rate[i])
        else
            irri_dem_gross = 0.0
        end
        model.variables.demand_gross[i] = irri_dem_gross
    end
    return nothing
end
update_demand_gross!(model::NoIrrigationNonPaddy, soil::SbmSoilModel) = nothing

"Struct to store paddy irrigation model variables"
@with_kw struct PaddyVariables
    n::Int
    demand_gross::Vector{Float64} = fill(MISSING_VALUE, n) # irrigation gross demand [mm Δt⁻¹]
    h::Vector{Float64} = zeros(n)                          # actual water depth in rice field [mm]
    evaporation::Vector{Float64} = zeros(n)                # evaporation rate [mm Δt⁻¹]
end

"Struct to store paddy irrigation model parameters"
@with_kw struct PaddyParameters
    irrigation_efficiency::Vector{Float64}        # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{Float64}      # maximum irrigation rate [mm Δt⁻¹]
    irrigation_areas::Vector{Bool}          # irrigation areas [-]
    irrigation_trigger::Vector{Bool}        # irrigation on or off [-]
    h_min::Vector{Float64}                        # minimum required water depth in the irrigated rice field [mm]
    h_opt::Vector{Float64}                        # optimal water depth in the irrigated rice fields [mm]
    h_max::Vector{Float64}                        # water depth when rice field starts spilling water (overflow) [mm]
end

"Paddy (flooded rice) irrigation model"
@with_kw struct Paddy <: AbstractIrrigationModel
    parameters::PaddyParameters
    variables::PaddyVariables
end

"Initialize paddy irrigation model"
function Paddy(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    h_min = ncread(
        dataset,
        config,
        "irrigated_paddy__min_depth";
        sel = indices,
        defaults = 20.0,
        type = Float64,
    )
    h_opt = ncread(
        dataset,
        config,
        "irrigated_paddy__optimal_depth";
        sel = indices,
        defaults = 50.0,
        type = Float64,
    )
    h_max = ncread(
        dataset,
        config,
        "irrigated_paddy__max_depth";
        sel = indices,
        defaults = 80.0,
        type = Float64,
    )
    efficiency = ncread(
        dataset,
        config,
        "irrigated_paddy__irrigation_efficiency";
        sel = indices,
        defaults = 1.0,
        type = Float64,
    )
    areas = ncread(
        dataset,
        config,
        "irrigated_paddy_area__count";
        optional = false,
        sel = indices,
        type = Bool,
    )
    irrigation_trigger = ncread(
        dataset,
        config,
        "irrigated_paddy__irrigation_trigger_flag";
        optional = false,
        sel = indices,
        type = Bool,
    )
    max_irri_rate =
        ncread(
            dataset,
            config,
            "irrigated_paddy__max_irrigation_rate";
            sel = indices,
            defaults = 25.0,
            type = Float64,
        ) .* (dt / BASETIMESTEP)
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
    variables = PaddyVariables(; n)
    paddy = Paddy(; parameters, variables)
    return paddy
end

# wrapper methods
get_water_depth(model::Paddy) = model.variables.h
get_water_depth(model::NoIrrigationPaddy) = Zeros(model.n)
get_demand_gross(model::Paddy) = model.variables.demand_gross
get_demand_gross(model::NoIrrigationPaddy) = Zeros(model.n)

"""
    evaporation!(model::Paddy, potential_evaporation)

Update `evaporation` and the water depth `h` of the paddy irrigation model for a single
timestep.
"""
function evaporation!(model::Paddy, potential_evaporation)
    for i in eachindex(potential_evaporation)
        if model.parameters.irrigation_areas[i]
            evaporation = min(model.variables.h[i], potential_evaporation[i])
            model.variables.h[i] -= evaporation
            model.variables.evaporation[i] = evaporation
        end
    end
    return nothing
end
evaporation!(model::NoIrrigationPaddy, potential_evaporation) = nothing

# wrapper methods
get_evaporation(model::NoIrrigationPaddy) = Zeros(model.n)
get_evaporation(model::Paddy) = model.variables.evaporation

"""
    update_runoff!(model::Paddy, runoff)

Update `runoff` based on the water depth `h_max` (paddy field starts spilling), and update
the water depth `h` of the paddy irrigation model for a single timestep.
"""
function update_runoff!(model::Paddy, runoff)
    for i in eachindex(model.parameters.irrigation_areas)
        if model.parameters.irrigation_areas[i]
            paddy_runoff = max(runoff[i] - model.parameters.h_max[i], 0.0)
            model.variables.h[i] = runoff[i] - paddy_runoff
            runoff[i] = paddy_runoff
        end
    end
    return nothing
end
update_runoff!(model::NoIrrigationPaddy, runoff) = nothing

"""
    update_demand_gross!(model::Paddy)

Update gross water demand `demand_gross` of the paddy irrigation model for a single
timestep.

The gross water demand is based on irrigation that is applied when the `irrigation_trigger`
is `true` (`on`) and when the paddy water depth `h` reaches below the minimum water depth
`h_min`. Irrigation is the amount required to reach the optimal paddy water depth `h_opt`,
taking into account limited irrigation efficiency and limited by a maximum irrigation rate.
"""
function update_demand_gross!(model::Paddy)
    (;
        irrigation_areas,
        irrigation_trigger,
        irrigation_efficiency,
        maximum_irrigation_rate,
        h_opt,
        h_min,
    ) = model.parameters
    (; h, demand_gross) = model.variables
    for i in eachindex(irrigation_areas)
        if irrigation_areas[i] && irrigation_trigger[i]
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
            irri_dem_gross = irr_depth_paddy / irrigation_efficiency[i]
            # limit irrigation demand to the maximum irrigation rate
            irri_dem_gross = min(irri_dem_gross, maximum_irrigation_rate[i])
        else
            irri_dem_gross = 0.0
        end
        demand_gross[i] = irri_dem_gross
    end
    return nothing
end

update_demand_gross!(model::NoIrrigationPaddy) = nothing

"Struct to store water demand model variables"
@with_kw struct DemandVariables
    n::Int
    irri_demand_gross::Vector{Float64} = zeros(n)        # irrigation gross demand [mm Δt⁻¹]
    nonirri_demand_gross::Vector{Float64} = zeros(n)     # non-irrigation gross demand [mm Δt⁻¹]
    total_gross_demand::Vector{Float64} = zeros(n)       # total gross demand [mm Δt⁻¹]
    surfacewater_demand::Vector{Float64} = zeros(n)      # demand from surface water [mm Δt⁻¹]
    groundwater_demand::Vector{Float64} = zeros(n)       # demand from groundwater [mm Δt⁻¹]
end

"Water demand model"
@with_kw struct Demand{
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

@with_kw struct NoDemand <: AbstractDemandModel
    n::Int
    domestic::NoNonIrrigationDemand = NoNonIrrigationDemand(n)
    industry::NoNonIrrigationDemand = NoNonIrrigationDemand(n)
    livestock::NoNonIrrigationDemand = NoNonIrrigationDemand(n)
    paddy::NoIrrigationPaddy = NoIrrigationPaddy(n)
    nonpaddy::NoIrrigationNonPaddy = NoIrrigationNonPaddy(n)
end

"Initialize water demand model"
function Demand(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    n = length(indices)
    demand(name; constr = NonIrrigationDemand, constr_triv = NoNonIrrigationDemand) =
        if getfield(config.model.water_demand, Symbol("$(name)__flag"))::Bool
            if constr == NonIrrigationDemand
                constr(dataset, config, indices, dt, name)
            else
                constr(dataset, config, indices, dt)
            end
        else
            constr_triv(n)
        end

    domestic = demand("domestic")
    industry = demand("industry")
    livestock = demand("livestock")
    paddy = demand("paddy"; constr = Paddy, constr_triv = NoIrrigationPaddy)
    nonpaddy = demand("nonpaddy"; constr = NonPaddy, constr_triv = NoIrrigationNonPaddy)

    variables = DemandVariables(; n)
    demand = Demand(; domestic, industry, livestock, paddy, nonpaddy, variables)
    return demand
end

"Struct to store river allocation model variables"
@with_kw struct AllocationRiverVariables
    n::Int
    act_surfacewater_abst::Vector{Float64} = zeros(n)        # actual surface water abstraction [mm Δt⁻¹]
    act_surfacewater_abst_vol::Vector{Float64} = zeros(n)    # actual surface water abstraction [m³ Δt⁻¹]
    available_surfacewater::Vector{Float64} = zeros(n)       # available surface water [m³]
    nonirri_returnflow::Vector{Float64} = zeros(n)           # return flow from non irrigation [mm Δt⁻¹]
end

"River allocation model"
@with_kw struct AllocationRiver <: AbstractAllocationModel
    n::Int
    variables::AllocationRiverVariables = AllocationRiverVariables(; n)
end

get_nonirrigation_returnflow(model::AllocationRiver) = model.variables.nonirri_returnflow
get_nonirrigation_returnflow(model::NoAllocationRiver) = Zeros(model.n)

"Struct to store land allocation allocation model parameters"
@with_kw struct AllocationLandParameters
    frac_sw_used::Vector{Float64}     # fraction surface water used [-]
    areas::Vector{Int}          # allocation areas [-]
end

"Struct to store land allocation model variables"
@with_kw struct AllocationLandVariables
    n::Int
    surfacewater_alloc::Vector{Float64} = zeros(n)           # allocation from surface water [mm Δt⁻¹]
    act_groundwater_abst::Vector{Float64} = zeros(n)         # actual groundwater abstraction [mm Δt⁻¹]
    act_groundwater_abst_vol::Vector{Float64} = zeros(n)     # actual groundwater abstraction [m³ Δt⁻¹]
    available_groundwater::Vector{Float64} = zeros(n)        # available groundwater [m³]
    groundwater_alloc::Vector{Float64} = zeros(n)            # allocation from groundwater [mm Δt⁻¹]
    irri_alloc::Vector{Float64} = zeros(n)                   # allocated water for irrigation [mm Δt⁻¹]
    nonirri_alloc::Vector{Float64} = zeros(n)                # allocated water for non-irrigation [mm Δt⁻¹]
    total_alloc::Vector{Float64} = zeros(n)                  # total allocated water [mm Δt⁻¹]
    nonirri_returnflow::Vector{Float64} = zeros(n)           # return flow from non irrigation [mm Δt⁻¹]
end

"Land allocation model"
@with_kw struct AllocationLand <: AbstractAllocationModel
    n::Int
    parameters::AllocationLandParameters
    variables::AllocationLandVariables = AllocationLandVariables(; n)
end

"Initialize water allocation for the land domain"
function AllocationLand(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    frac_sw_used = ncread(
        dataset,
        config,
        "land_surface_water__withdrawal_fraction";
        sel = indices,
        defaults = 1,
        type = Float64,
    )
    areas = ncread(
        dataset,
        config,
        "land_water_allocation_area__count";
        sel = indices,
        defaults = 1,
        type = Int,
    )

    n = length(indices)

    parameters = AllocationLandParameters(; areas, frac_sw_used)
    allocation = AllocationLand(; n, parameters)
    return allocation
end

# wrapper methods
get_irrigation_allocated(model::AllocationLand) = model.variables.irri_alloc
get_irrigation_allocated(model::NoAllocationLand) = Zeros(model.n)
get_nonirrigation_returnflow(model::AllocationLand) = model.variables.nonirri_returnflow
get_nonirrigation_returnflow(model::NoAllocationLand) = Zeros(model.n)
get_groundwater_abstraction_flux(model::AllocationLand) =
    model.variables.act_groundwater_abst
get_groundwater_abstraction_flux(model::NoAllocationLand) = Zeros(model.n)

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
function return_flow_fraction!(model::NonIrrigationDemand)
    (; returnflow_fraction) = model.variables
    (; demand_gross, demand_net) = model.demand
    @. returnflow_fraction = return_flow_fraction(demand_gross, demand_net)
    return nothing
end

# return zero (gross water demand) if non-irrigation water demand sector is not defined
return_flow_fraction!(model::NoNonIrrigationDemand) = nothing

"""
Update water allocation for river and land domains based on local surface water (river)
availability.
"""
function surface_water_allocation_local!(
    model::AllocationLand,
    demand::Demand,
    river::AbstractRiverFlowModel,
    domain::DomainLand,
    dt::Float64,
)
    (; surfacewater_alloc) = model.variables
    (; surfacewater_demand) = demand.variables
    (; act_surfacewater_abst_vol, act_surfacewater_abst, available_surfacewater) =
        river.allocation.variables
    (; external_inflow) = river.boundary_conditions
    (; storage) = river.variables

    (; area) = domain.parameters
    indices_river = domain.network.river_inds_excl_reservoir

    # maps from the land domain to the internal river domain (linear index), excluding reservoirs
    for i in eachindex(surfacewater_demand)
        if indices_river[i] > 0.0
            # the available volume is limited by a fixed scaling factor of 0.8 to prevent
            # rivers completely drying out. check for abstraction through negative external
            # inflow first and adjust available volume.
            if external_inflow[indices_river[i]] < 0.0
                available_volume = storage[indices_river[i]] * 0.80
                max_river_abstraction =
                    min(-external_inflow[indices_river[i]] * dt, available_volume)
                available_volume = max(available_volume - max_river_abstraction, 0.0)
            else
                available_volume = storage[indices_river[i]] * 0.80
            end
            # satisfy surface water demand with available local river volume
            surfacewater_demand_vol = surfacewater_demand[i] * 0.001 * area[i]
            abstraction_vol = min(surfacewater_demand_vol, available_volume)
            act_surfacewater_abst_vol[indices_river[i]] = abstraction_vol
            # remaining available surface water and demand
            available_surfacewater[indices_river[i]] =
                max(available_volume - abstraction_vol, 0.0)
            abstraction = (abstraction_vol / area[i]) * 1000.0
            surfacewater_demand[i] = max(surfacewater_demand[i] - abstraction, 0.0)
            # update actual abstraction from river and surface water allocation (land cell)
            act_surfacewater_abst[indices_river[i]] = abstraction
            surfacewater_alloc[i] = abstraction
        end
    end
    return nothing
end

"""
Update water allocation for river and land domains based on surface water (river)
availability for allocation areas.
"""
function surface_water_allocation_area!(
    model::AllocationLand,
    demand::Demand,
    river::AbstractRiverFlowModel,
    domain::Domain,
    dt::Float64,
)
    inds_river = domain.river.network.allocation_area_indices
    inds_land = domain.land.network.allocation_area_indices
    inds_reservoir = domain.river.network.reservoir_indices
    (; area) = domain.land.parameters

    (; available_surfacewater, act_surfacewater_abst_vol, act_surfacewater_abst) =
        river.allocation.variables
    (; surfacewater_alloc) = model.variables
    (; surfacewater_demand) = demand.variables
    (; reservoir) = river.boundary_conditions

    # loop over allocation areas
    for i in eachindex(inds_river)
        # surface water demand (allocation area)
        sw_demand_vol = 0.0
        for j in inds_land[i]
            sw_demand_vol += surfacewater_demand[j] * 0.001 * area[j]
        end
        # surface water availability (allocation area)
        sw_available = 0.0
        for j in inds_river[i]
            if inds_reservoir[j] > 0
                # for reservoir locations use reservoir storage, check for abstraction
                # through external negative inflow first and adjust available volume.
                k = inds_reservoir[j]
                external_inflow = reservoir.boundary_conditions.external_inflow[k]
                if external_inflow < 0.0
                    available_volume = reservoir.variables.storage[k] * 0.98
                    max_res_abstraction = min(-external_inflow * dt, available_volume)
                    available_volume = max(available_volume - max_res_abstraction, 0.0)
                else
                    available_volume = reservoir.variables.storage[k] * 0.98
                end
                available_surfacewater[j] = available_volume
                sw_available += available_volume
            else
                # river volume
                sw_available += available_surfacewater[j]
            end
        end
        # total actual surface water abstraction [m3] in an allocation area, minimum of
        # available surface water and demand in an allocation area.
        sw_abstraction = min(sw_available, sw_demand_vol)

        # fraction of available surface water that can be abstracted at allocation area
        # level
        frac_abstract_sw = bounded_divide(sw_abstraction, sw_available)
        # fraction of water demand that can be satisfied by available surface water at
        # allocation area level.
        frac_allocate_sw = bounded_divide(sw_abstraction, sw_demand_vol)

        # water abstracted from surface water at each river cell (including reservoir
        # locations).
        for j in inds_river[i]
            act_surfacewater_abst_vol[j] += frac_abstract_sw * available_surfacewater[j]
            act_surfacewater_abst[j] =
                (act_surfacewater_abst_vol[j] / domain.river.parameters.cell_area[j]) *
                1000.0
        end

        # water allocated to each land cell.
        for j in inds_land[i]
            surfacewater_alloc[j] += frac_allocate_sw * surfacewater_demand[j]
        end
    end
    return nothing
end

"Update water allocation for land domain based on local groundwater availability."
function groundwater_allocation_local!(
    model::AllocationLand,
    demand::Demand,
    groundwater_storage::Vector{Float64},
    parameters::LandParameters,
)
    (;
        surfacewater_alloc,
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = model.variables
    (; groundwater_demand, total_gross_demand) = demand.variables
    (; area, reservoir_coverage) = parameters

    for i in eachindex(groundwater_demand)
        # groundwater demand based on allocation from surface water.
        groundwater_demand[i] = max(total_gross_demand[i] - surfacewater_alloc[i], 0.0)
        # excluding reservoirs
        if !reservoir_coverage[i]
            # satisfy groundwater demand with available local groundwater volume
            groundwater_demand_vol = groundwater_demand[i] * 0.001 * area[i]
            available_volume = groundwater_storage[i] * 0.75 # limit available groundwater volume
            abstraction_vol = min(groundwater_demand_vol, available_volume)
            act_groundwater_abst_vol[i] = abstraction_vol
            # remaining available groundwater and demand
            available_groundwater[i] = max(available_volume - abstraction_vol, 0.0)
            abstraction = (abstraction_vol / area[i]) * 1000.0
            groundwater_demand[i] = max(groundwater_demand[i] - abstraction, 0.0)
            # update actual abstraction from groundwater and groundwater allocation (land cell)
            act_groundwater_abst[i] = abstraction
            groundwater_alloc[i] = abstraction
        end
    end
    return nothing
end

"""
Update water allocation for land domain based on groundwater availability for allocation
areas.

"""
function groundwater_allocation_area!(model::AllocationLand, demand::Demand, domain::Domain)
    inds_river = domain.river.network.allocation_area_indices
    inds_land = domain.land.network.allocation_area_indices
    (;
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = model.variables

    (; groundwater_demand) = demand.variables
    (; area) = domain.land.parameters

    # loop over allocation areas
    for i in eachindex(inds_river)
        # groundwater demand and availability (allocation area)
        gw_demand_vol = 0.0
        gw_available = 0.0
        for j in inds_land[i]
            gw_demand_vol += groundwater_demand[j] * 0.001 * area[j]
            gw_available += available_groundwater[j]
        end
        # total actual groundwater abstraction [m3] in an allocation area, minimum of
        # available  groundwater and demand in an allocation area.
        gw_abstraction = min(gw_available, gw_demand_vol)

        # fraction of available groundwater that can be abstracted at allocation area level
        frac_abstract_gw = bounded_divide(gw_abstraction, gw_available)
        # fraction of water demand that can be satisfied by available groundwater at
        # allocation area level.
        frac_allocate_gw = bounded_divide(gw_abstraction, gw_demand_vol)

        # water abstracted from groundwater and allocated.
        for j in inds_land[i]
            act_groundwater_abst_vol[j] += frac_abstract_gw * available_groundwater[j]
            act_groundwater_abst[j] = 1000.0 * (act_groundwater_abst_vol[j] / area[j])
            groundwater_alloc[j] += frac_allocate_gw * groundwater_demand[j]
        end
    end
    return nothing
end

"Return and update non-irrigation sector (domestic, livestock, industry) return flow"
function return_flow(
    model::NonIrrigationDemand,
    nonirri_demand_gross::Vector{Float64},
    nonirri_alloc::Vector{Float64},
)
    for i in eachindex(model.variables.returnflow)
        frac = bounded_divide(model.demand.demand_gross[i], nonirri_demand_gross[i])
        allocate = frac * nonirri_alloc[i]
        model.variables.returnflow[i] = model.variables.returnflow_fraction[i] * allocate
    end
    return model.variables.returnflow
end

# return zero (return flow) if non-irrigation sector is not defined
return_flow(
    model::NoNonIrrigationDemand,
    nonirri_demand_gross::Vector{Float64},
    nonirri_alloc::Vector{Float64},
) = 0.0

# wrapper methods
groundwater_storage(model::LateralSSF) = model.variables.storage
groundwater_storage(model) = model.aquifer.variables.storage

"""
    update_water_allocation!(
    model::AllocationLand,
    demand::Demand,
    routing::Routing,
    domain::Domain,
    dt::Float64,
)

Update water allocation for the land domain `AllocationLand` and water allocation for the
river domain (part of `routing`) based on the water `demand` model for a single timestep.
First, surface water abstraction is computed to satisfy local water demand (non-irrigation
and irrigation), and then updated (including reservoirs) to satisfy the remaining water
demand for allocation areas. Then groundwater abstraction is computed to satisfy the
remaining local water demand, and then updated to satisfy the remaining water demand for
allocation areas. Finally, non-irrigation return flows are updated.
"""
function update_water_allocation!(
    model::AllocationLand,
    demand::Demand,
    routing::Routing,
    domain::Domain,
    dt::Float64,
)
    river = routing.river_flow
    inds_river = domain.land.network.river_inds_excl_reservoir
    inds_reservoir = domain.reservoir.network.river_indices
    (;
        groundwater_alloc,
        surfacewater_alloc,
        act_groundwater_abst,
        act_groundwater_abst_vol,
        total_alloc,
        irri_alloc,
        nonirri_alloc,
        nonirri_returnflow,
    ) = model.variables

    (; surfacewater_demand, nonirri_demand_gross, irri_demand_gross, total_gross_demand) =
        demand.variables

    (; frac_sw_used) = model.parameters
    (; act_surfacewater_abst, act_surfacewater_abst_vol) = river.allocation.variables
    (; abstraction, reservoir) = river.boundary_conditions

    surfacewater_alloc .= 0.0
    act_surfacewater_abst .= 0.0
    act_surfacewater_abst_vol .= 0.0
    # total surface water demand for each land cell
    @. surfacewater_demand =
        frac_sw_used * nonirri_demand_gross + frac_sw_used * irri_demand_gross

    # local surface water demand and allocation (river, excluding reservoirs)
    surface_water_allocation_local!(model, demand, river, domain.land, dt)
    # surface water demand and allocation for areas
    surface_water_allocation_area!(model, demand, river, domain, dt)

    @. abstraction = act_surfacewater_abst_vol / dt

    # for reservoir locations set river abstraction at zero and abstract volume
    # from reservoir, including an update of waterlevel
    if !isnothing(reservoir)
        @. abstraction[inds_reservoir] = 0.0
        @. reservoir.variables.storage -= act_surfacewater_abst_vol[inds_reservoir]
        @. reservoir.variables.waterlevel = waterlevel(
            reservoir.parameters.storfunc,
            reservoir.parameters.area,
            reservoir.variables.storage,
            reservoir.parameters.sh,
        )
    end

    groundwater_alloc .= 0.0
    act_groundwater_abst_vol .= 0.0
    act_groundwater_abst .= 0.0
    # local groundwater demand and allocation
    groundwater_allocation_local!(
        model,
        demand,
        groundwater_storage(routing.subsurface_flow),
        domain.land.parameters,
    )
    # groundwater demand and allocation for areas
    groundwater_allocation_area!(model, demand, domain)

    # irrigation allocation
    for i in eachindex(total_alloc)
        total_alloc[i] = groundwater_alloc[i] + surfacewater_alloc[i]
        frac_irri = bounded_divide(irri_demand_gross[i], total_gross_demand[i])
        irri_alloc[i] = frac_irri * total_alloc[i]
        nonirri_alloc[i] = total_alloc[i] - irri_alloc[i]
    end

    # non-irrigation return flows
    returnflow_livestock =
        return_flow(demand.livestock, nonirri_demand_gross, nonirri_alloc)
    returnflow_domestic = return_flow(demand.domestic, nonirri_demand_gross, nonirri_alloc)
    returnflow_industry = return_flow(demand.industry, nonirri_demand_gross, nonirri_alloc)

    @. nonirri_returnflow = returnflow_livestock + returnflow_domestic + returnflow_industry

    for i in eachindex(nonirri_returnflow)
        if inds_river[i] > 0.0
            k = inds_river[i]
            river.allocation.variables.nonirri_returnflow[k] = nonirri_returnflow[i]
            nonirri_returnflow[i] = 0.0
        else
            nonirri_returnflow[i] = nonirri_returnflow[i]
        end
    end
end
update_water_allocation!(
    model::NoAllocationLand,
    demand::NoDemand,
    routing::Routing,
    domain::Domain,
    dt::Float64,
) = nothing

"""
    update_demand_gross!(model::Demand)

Update total irrigation gross water demand `irri_demand_gross`, total non-irrigation gross
water demand `nonirri_demand_gross` and total gross water demand `total_gross_demand`.
"""
function update_demand_gross!(model::Demand)
    (; nonpaddy, paddy, domestic, industry, livestock) = model
    (; irri_demand_gross, nonirri_demand_gross, total_gross_demand) = model.variables
    # get gross water demands
    industry_dem = get_demand_gross(industry)
    domestic_dem = get_demand_gross(domestic)
    livestock_dem = get_demand_gross(livestock)
    nonpaddy_dem_gross = get_demand_gross(nonpaddy)
    paddy_dem_gross = get_demand_gross(paddy)
    # update gross water demands
    @. irri_demand_gross = nonpaddy_dem_gross + paddy_dem_gross
    @. nonirri_demand_gross = industry_dem + domestic_dem + livestock_dem
    @. total_gross_demand =
        nonpaddy_dem_gross + paddy_dem_gross + industry_dem + domestic_dem + livestock_dem

    return nothing
end

update_demand_gross!(model::NoDemand) = nothing

"""
    update_water_demand!(model::Demand, soil::SbmSoilModel)

Update the return flow fraction `returnflow_fraction` of `industry`, `domestic` and
`livestock`, gross water demand `demand_gross` of `paddy` and `nonpaddy` models, and the
total gross water demand, total irrigation gross water demand and total non-irrigation gross
water demand as part of the water `demand` model.
"""
function update_water_demand!(model::Demand, soil::SbmSoilModel)
    (; nonpaddy, paddy, domestic, industry, livestock) = model

    return_flow_fraction!(industry)
    return_flow_fraction!(domestic)
    return_flow_fraction!(livestock)

    update_demand_gross!(nonpaddy, soil)
    update_demand_gross!(paddy)
    update_demand_gross!(model)

    return nothing
end
update_water_demand!(model::NoDemand, soil::SbmSoilModel) = nothing
