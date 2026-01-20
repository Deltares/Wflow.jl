abstract type AbstractIrrigationModel end
abstract type AbstractIrrigationDemandModel <: AbstractDemandModel end

struct NoIrrigationPaddy <: AbstractIrrigationModel
    n::Int
end
struct NoIrrigationNonPaddy <: AbstractIrrigationModel
    n::Int
end
struct NoNonIrrigationDemand <: AbstractIrrigationDemandModel
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
    # return flow [mm dt⁻¹ => m s⁻¹]
    returnflow::Vector{Float64}
    # return flow fraction [-]
    returnflow_fraction::Vector{Float64}
end

"Struct to store prescribed water demand variables"
@with_kw struct PrescibedDemand
    # gross water demand [mm dt⁻¹ => m s⁻¹]
    demand_gross::Vector{Float64}
    # net water demand [mm dt⁻¹ => m s⁻¹]
    demand_net::Vector{Float64}
end

"Non-irrigation water demand model"
@with_kw struct NonIrrigationDemand <: AbstractIrrigationDemandModel
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
    sector::AbstractString,
)
    demand_gross = ncread(
        dataset,
        config,
        "$(sector)__gross_water_demand_volume_flux",
        LandHydrologySBM;
        sel = indices,
        defaults = 0.0,
        type = Float64,
    )
    demand_net = ncread(
        dataset,
        config,
        "$(sector)__net_water_demand_volume_flux",
        LandHydrologySBM;
        sel = indices,
        defaults = 0.0,
        type = Float64,
    )
    n = length(indices)
    returnflow_f = return_flow_fraction.(demand_gross, demand_net)

    demand = PrescibedDemand(; demand_gross, demand_net)
    vars = NonIrrigationDemandVariables(;
        returnflow_fraction = returnflow_f,
        returnflow = zeros(n),
    )
    non_irrigation_demand = NonIrrigationDemand(; demand, variables = vars)

    return non_irrigation_demand
end

"Struct to store non-paddy irrigation model variables"
@with_kw struct NonPaddyVariables
    n::Int
    demand_gross::Vector{Float64} = fill(MISSING_VALUE, n)   # irrigation gross demand [mm dt⁻¹ => m s⁻¹]
end

"Struct to store non-paddy irrigation model parameters"
@with_kw struct NonPaddyParameters
    irrigation_efficiency::Vector{Float64}        # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{Float64}      # maximum irrigation rate [mm dt⁻¹ => m s⁻¹]
    irrigation_areas::Vector{Bool}                # irrigation areas [-]
    irrigation_trigger::Vector{Bool}              # irrigation on or off [-]
end

"Non-paddy (other crops than flooded rice) irrigation model"
@with_kw struct NonPaddy <: AbstractIrrigationModel
    n::Int
    parameters::NonPaddyParameters
    variables::NonPaddyVariables = NonPaddyVariables(; n)
end

"Initialize non-paddy irrigation model"
function NonPaddy(dataset::NCDataset, config::Config, indices::Vector{CartesianIndex{2}})
    efficiency = ncread(
        dataset,
        config,
        "irrigated_non_paddy__irrigation_efficiency",
        LandHydrologySBM;
        sel = indices,
        defaults = 1.0,
        type = Float64,
    )
    areas = ncread(
        dataset,
        config,
        "irrigated_non_paddy_area__count",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Int,
    )
    irrigation_trigger = ncread(
        dataset,
        config,
        "irrigated_non_paddy__irrigation_trigger_flag",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Bool,
    )
    max_irri_rate = ncread(
        dataset,
        config,
        "irrigated_non_paddy__max_irrigation_rate",
        LandHydrologySBM;
        sel = indices,
        defaults = 25.0,
        type = Float64,
    )

    parameters = NonPaddyParameters(;
        maximum_irrigation_rate = max_irri_rate,
        irrigation_efficiency = efficiency,
        irrigation_areas = areas,
        irrigation_trigger,
    )
    n = length(indices)
    nonpaddy = NonPaddy(; n, parameters)

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
    (; irrigation_areas, irrigation_trigger, maximum_irrigation_rate) = model.parameters
    (; demand_gross) = model.variables
    (; n_unsatlayers) = soil.variables

    for i in eachindex(irrigation_areas)
        if irrigation_areas[i] && irrigation_trigger[i]
            # [m]
            irri_dem_gross_depth = 0.0
            for k in 1:n_unsatlayers[i]
                depletion, readily_available_water = water_demand_root_zone(soil, i, k)

                # check if maximum irrigation rate has been applied at the previous time step.
                max_irri_rate_applied =
                    model.variables.demand_gross[i] == maximum_irrigation_rate[i]
                if depletion >= readily_available_water # start irrigation
                    # [m] += [m]
                    irri_dem_gross += depletion
                    # add depletion to irrigation gross demand when the maximum irrigation rate has been
                    # applied at the previous time step (to get volumetric water content at field capacity)
                elseif depletion > 0.0 && max_irri_rate_applied # continue irrigation
                    # [m] += [m]
                    irri_dem_gross_depth += depletion
                end
            end
            demand_gross[i] = compute_demand_gross(model, soil, irri_dem_gross, i)
        else
            demand_gross[i] = 0.0
        end
    end
    return nothing
end

update_demand_gross!(model::NoIrrigationNonPaddy, soil::SbmSoilModel) = nothing

"Compute water demand only for root zone through root fraction per layer"
function water_demand_root_zone(soil::SbmSoilModel, i::Int, k::Int)
    (; sumlayers, hb, theta_s, theta_r, c) = soil.parameters
    (; ustorelayerthickness, ustorelayerdepth, h3) = soil.variables

    rootingdepth = get_rootingdepth(soil)

    # [-] = clamp(([m] - [m]) / [m], [-], [-])
    rootfrac =
        clamp((rootingdepth[i] - sumlayers[i][k]) / ustorelayerthickness[i][k], 0.0, 1.0)
    # vwc_f and vwc_h3 can be precalculated.
    # [-]
    vwc_fc = vwc_brooks_corey(-1.0, hb[i], theta_s[i], theta_r[i], c[i][k])
    # [-]
    vwc_h3 = vwc_brooks_corey(h3[i], hb[i], theta_s[i], theta_r[i], c[i][k])
    # [m] = ([-] * [m]) - ([m] + [-] * [m])
    depletion =
        (vwc_fc * ustorelayerthickness[i][k]) -
        (ustorelayerdepth[i][k] + theta_r[i] * ustorelayerthickness[i][k])
    # [m] = [m] * [-]
    depletion *= rootfrac
    # [m] = [-] * [m]
    readily_available_water = (vwc_fc - vwc_h3) * ustorelayerthickness[i][k]
    # [m] = [m] * [-]
    readily_available_water *= rootfrac

    return depletion, readily_available_water
end

function compute_demand_gross(
    model::NonPaddy,
    soil::SbmSoilModel,
    irri_dem_gross::Float64,
    i::Int,
)
    (; pathfrac, infiltcapsoil) = soil.parameters
    (; f_infiltration_reduction) = soil.variables
    (; irrigation_efficiency, maximum_irrigation_rate) = model.parameters

    # [m s⁻¹] = [-] * [-] * [m s⁻¹]
    infiltration_capacity =
        f_infiltration_reduction[i] * (1.0 - pathfrac[i]) * infiltcapsoil[i]
    # [m s⁻¹] = min([m s⁻¹], [m s⁻¹])
    irri_dem_gross = min(irri_dem_gross, infiltration_capacity)
    # [m s⁻¹] = [m s⁻¹] / [-]
    irri_dem_gross /= irrigation_efficiency[i]
    # limit irrigation demand to the maximum irrigation rate
    irri_dem_gross = min(irri_dem_gross, maximum_irrigation_rate[i])

    return irri_dem_gross
end

"Struct to store paddy irrigation model variables"
@with_kw struct PaddyVariables
    n::Int
    demand_gross::Vector{Float64} = fill(MISSING_VALUE, n) # irrigation gross demand [mm dt⁻¹ => m s⁻¹]
    h::Vector{Float64} = zeros(n)                          # actual water depth in rice field [mm => m]
    evaporation::Vector{Float64} = zeros(n)                # evaporation rate [mm dt⁻¹ => m s⁻¹]
end

"Struct to store paddy irrigation model parameters"
@with_kw struct PaddyParameters
    irrigation_efficiency::Vector{Float64}        # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{Float64}      # maximum irrigation rate [mm dt⁻¹ => m s⁻¹]
    irrigation_areas::Vector{Bool}                # irrigation areas [-]
    irrigation_trigger::Vector{Bool}              # irrigation on or off [-]
    h_min::Vector{Float64}                        # minimum required water depth in the irrigated rice field [mm => m]
    h_opt::Vector{Float64}                        # optimal water depth in the irrigated rice fields  => m
    h_max::Vector{Float64}                        # water depth when rice field starts spilling water (overflow)  => m
end

"Paddy (flooded rice) irrigation model"
@with_kw struct Paddy <: AbstractIrrigationModel
    n::Int
    parameters::PaddyParameters
    variables::PaddyVariables = PaddyVariables(; n)
end

"Initialize paddy irrigation model"
function Paddy(dataset::NCDataset, config::Config, indices::Vector{CartesianIndex{2}})
    h_min = ncread(
        dataset,
        config,
        "irrigated_paddy__min_depth",
        LandHydrologySBM;
        sel = indices,
        defaults = 20.0,
        type = Float64,
    )
    h_opt = ncread(
        dataset,
        config,
        "irrigated_paddy__optimal_depth",
        LandHydrologySBM;
        sel = indices,
        defaults = 50.0,
        type = Float64,
    )
    h_max = ncread(
        dataset,
        config,
        "irrigated_paddy__max_depth",
        LandHydrologySBM;
        sel = indices,
        defaults = 80.0,
        type = Float64,
    )
    efficiency = ncread(
        dataset,
        config,
        "irrigated_paddy__irrigation_efficiency",
        LandHydrologySBM;
        sel = indices,
        defaults = 1.0,
        type = Float64,
    )
    areas = ncread(
        dataset,
        config,
        "irrigated_paddy_area__count",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Bool,
    )
    irrigation_trigger = ncread(
        dataset,
        config,
        "irrigated_paddy__irrigation_trigger_flag",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Bool,
    )
    max_irri_rate = ncread(
        dataset,
        config,
        "irrigated_paddy__max_irrigation_rate",
        LandHydrologySBM;
        sel = indices,
        defaults = 25.0,
        type = Float64,
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
    paddy = Paddy(; n, parameters)
    return paddy
end

# wrapper methods
get_water_depth(model::Paddy) = model.variables.h
get_water_depth(model::NoIrrigationPaddy) = Zeros(model.n)
get_demand_gross(model::Paddy) = model.variables.demand_gross
get_demand_gross(model::NoIrrigationPaddy) = Zeros(model.n)

"""
    evaporation!(model::Paddy, potential_evaporation, dt::Number)

Update `evaporation` and the water depth `h` of the paddy irrigation model for a single
timestep.
"""
function evaporation!(model::Paddy, potential_evaporation, dt::Number)
    for i in eachindex(potential_evaporation)
        if model.parameters.irrigation_areas[i]
            # [m s⁻¹] = min([h] / [s], [m s⁻¹])
            evaporation = min(model.variables.h[i] / dt, potential_evaporation[i])
            model.variables.h[i] -= evaporation * dt
            model.variables.evaporation[i] = evaporation
        end
    end
    return nothing
end
evaporation!(model::NoIrrigationPaddy, potential_evaporation, dt::Number) = nothing

# wrapper methods
get_evaporation(model::NoIrrigationPaddy) = Zeros(model.n)
get_evaporation(model::Paddy) = model.variables.evaporation

"""
    update_runoff!(model::Paddy, runoff, dt::Number)

Update `runoff` based on the water depth `h_max` (paddy field starts spilling), and update
the water depth `h` of the paddy irrigation model for a single timestep.
"""
function update_runoff!(model::Paddy, runoff, dt::Number)
    for i in eachindex(model.parameters.irrigation_areas)
        if model.parameters.irrigation_areas[i]
            # [m s⁻¹] = max([m s⁻¹] - [m] / [s], 0.0)
            paddy_runoff = max(runoff[i] - model.parameters.h_max[i] / dt, 0.0)
            # [m] = [m s⁻¹] * [s]
            model.variables.h[i] = (runoff[i] - paddy_runoff) * dt
            runoff[i] = paddy_runoff
        end
    end
    return nothing
end
update_runoff!(model::NoIrrigationPaddy, runoff, dt::Number) = nothing

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
    (; demand_gross) = model.variables
    (;
        irrigation_areas,
        irrigation_trigger,
        irrigation_efficiency,
        maximum_irrigation_rate,
    ) = model.parameters

    for i in eachindex(irrigation_areas)
        if irrigation_areas[i] && irrigation_trigger[i]
            # [m]
            irr_depth_paddy = compute_irrigation_depth(model, i)
            # [m] = [m] / [-]
            irri_dem_gross = irr_depth_paddy / irrigation_efficiency[i]
            # limit irrigation demand to the maximum irrigation rate
            # [m s⁻¹] = min([m] / [s], [m s⁻¹])
            irri_dem_gross = min(irri_dem_gross, maximum_irrigation_rate[i])
            demand_gross[i] = irri_dem_gross
        else
            demand_gross[i] = 0.0
        end
    end
end

function compute_irrigation_depth(model::Paddy, i::Int)
    (; maximum_irrigation_rate, h_min, h_opt) = model.parameters
    (; demand_gross, h) = model.variables

    # check if maximum irrigation rate has been applied at the previous time step.
    max_irri_rate_applied = demand_gross[i] == maximum_irrigation_rate[i]
    # start irrigation
    # [m]
    irr_depth_paddy = if h[i] < h_min[i]
        h_opt[i] - h[i]
    elseif h[i] < h_opt[i] && max_irri_rate_applied # continue irrigation
        h_opt[i] - h[i]
    else
        0.0
    end

    return irr_depth_paddy
end

update_demand_gross!(model::NoIrrigationPaddy) = nothing

"Struct to store water demand model variables"
@with_kw struct DemandVariables
    n::Int
    irri_demand_gross::Vector{Float64} = zeros(n)        # irrigation gross demand [mm dt⁻¹ => m s⁻¹]
    nonirri_demand_gross::Vector{Float64} = zeros(n)     # non-irrigation gross demand [mm dt⁻¹ => m s⁻¹]
    total_gross_demand::Vector{Float64} = zeros(n)       # total gross demand [mm dt⁻¹ => m s¹]
    surfacewater_demand::Vector{Float64} = zeros(n)      # demand from surface water [mm dt⁻¹ => m s⁻¹]
    groundwater_demand::Vector{Float64} = zeros(n)       # demand from groundwater [mm dt⁻¹ => m s⁻¹]
end

"Water demand model"
@with_kw struct Demand{
    D <: AbstractIrrigationDemandModel,
    I <: AbstractIrrigationDemandModel,
    L <: AbstractIrrigationDemandModel,
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
function Demand(dataset::NCDataset, config::Config, indices::Vector{CartesianIndex{2}})
    n = length(indices)
    demand(name; constr = NonIrrigationDemand, constr_triv = NoNonIrrigationDemand) =
        if getfield(config.model.water_demand, Symbol("$(name)__flag"))::Bool
            if constr == NonIrrigationDemand
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
    paddy = demand("paddy"; constr = Paddy, constr_triv = NoIrrigationPaddy)
    nonpaddy = demand("nonpaddy"; constr = NonPaddy, constr_triv = NoIrrigationNonPaddy)

    variables = DemandVariables(; n)
    return Demand(; domestic, industry, livestock, paddy, nonpaddy, variables)
end

"Struct to store river allocation model variables"
@with_kw struct AllocationRiverVariables
    n::Int
    act_surfacewater_abst::Vector{Float64} = zeros(n)       # actual surface water abstraction [mm dt⁻¹ => m s⁻¹]
    act_surfacewater_abst_vol::Vector{Float64} = zeros(n)   # actual surface water abstraction [m³ dt⁻¹ => m³ s⁻¹]
    available_surfacewater::Vector{Float64} = zeros(n)      # available surface water [m³]
    nonirri_returnflow::Vector{Float64} = zeros(n)          # return flow from non irrigation [mm dt⁻¹ => m s⁻¹]
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
    # allocation from surface water [mm dt⁻¹ => m s⁻¹]
    surfacewater_alloc::Vector{Float64} = zeros(n)
    # actual groundwater abstraction [mm dt⁻¹ => m s⁻¹]
    act_groundwater_abst::Vector{Float64} = zeros(n)
    # actual groundwater abstraction [m³ dt⁻¹ => m³ s⁻¹]
    act_groundwater_abst_vol::Vector{Float64} = zeros(n)
    # available groundwater [m³]
    available_groundwater::Vector{Float64} = zeros(n)
    # allocation from groundwater [mm dt⁻¹ => m s⁻¹]
    groundwater_alloc::Vector{Float64} = zeros(n)
    # allocated water for irrigation [mm dt⁻¹ => m s⁻¹]
    irri_alloc::Vector{Float64} = zeros(n)
    # allocated water for non-irrigation [mm dt⁻¹ => m s⁻¹]
    nonirri_alloc::Vector{Float64} = zeros(n)
    # total allocated water [mm dt⁻¹ => m s⁻¹]
    total_alloc::Vector{Float64} = zeros(n)
    # return flow from non irrigation [mm dt⁻¹ => m s⁻¹]
    nonirri_returnflow::Vector{Float64} = zeros(n)
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
        "land_surface_water__withdrawal_fraction",
        LandHydrologySBM;
        sel = indices,
        defaults = 1,
    )
    areas = ncread(
        dataset,
        config,
        "land_water_allocation_area__count",
        LandHydrologySBM;
        sel = indices,
        defaults = 1,
        type = Int,
    )

    n = length(indices)
    parameters = AllocationLandParameters(; areas, frac_sw_used)
    return AllocationLand(; n, parameters)
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
    demand_variables::DemandVariables,
    river::AbstractRiverFlowModel,
    domain::DomainLand,
    dt::Float64,
)
    (; surfacewater_alloc) = model.variables
    (; surfacewater_demand) = demand_variables
    (; act_surfacewater_abst_vol, act_surfacewater_abst, available_surfacewater) =
        river.allocation.variables
    (; external_inflow) = river.boundary_conditions
    (; storage) = river.variables
    (; area) = domain.parameters
    indices_river = domain.network.river_inds_excl_reservoir

    # maps from the land domain to the internal river domain (linear index), excluding reservoirs
    for (i, index_river) in enumerate(indices_river)
        if index_river > 0
            # the available volume is limited by a fixed scaling factor of 0.8 to prevent
            # rivers completely drying out. check for abstraction through negative external
            # inflow first and adjust available volume.
            if external_inflow[index_river] < 0.0
                # [m³] = [m³] * [-]
                available_volume = storage[index_river] * 0.80
                # [m³] = min([m³ s⁻¹] * [s], [m³])
                max_river_abstraction =
                    min(-external_inflow[index_river] * dt, available_volume)
                available_volume = max(available_volume - max_river_abstraction, 0.0)
            else
                # [m³] = [m³] * [-]
                available_volume = storage[index_river] * 0.80
            end
            # satisfy surface water demand with available local river volume
            # [m³ s⁻¹] = [m s⁻¹] * [m²]
            surfacewater_demand_vol = surfacewater_demand[i] * area[i]
            # [m³] = min([m³ s⁻¹] * [s], [m³])
            abstraction_vol = min(surfacewater_demand_vol * dt, available_volume)
            # [m³ s⁻¹] = [m³] / [s]
            act_surfacewater_abst_vol[index_river] = abstraction_vol / dt
            # remaining available surface water and demand [m³]
            available_surfacewater[index_river] =
                max(available_volume - abstraction_vol, 0.0)
            # [m s⁻¹] = [m³] / ([m²] * [s])
            abstraction = abstraction_vol / (area[i] * dt)
            # [m s⁻¹] = max([m s⁻¹] - [m s⁻¹], [m s⁻¹])
            surfacewater_demand[i] = max(surfacewater_demand[i] - abstraction, 0.0)
            # [m s⁻¹]
            # update actual abstraction from river and surface water allocation (land cell)
            act_surfacewater_abst[index_river] = abstraction
            # [m s⁻¹]
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
    demand_variables::DemandVariables,
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
    (; surfacewater_demand) = demand_variables
    (; reservoir) = river.boundary_conditions

    for i in eachindex(inds_river)
        # surface water_demand (allocation area)
        # [m³ s⁻¹] = ∑ [m s⁻¹] * [m²]
        sw_demand_vol = mapreduce(j -> surfacewater_demand[j] * area[j], +, inds_land[i])

        # [m³]
        sw_available = available_surface_water!(
            available_surfacewater,
            reservoir,
            inds_river[i],
            inds_reservoir,
            dt,
        )

        # total actual surface water abstraction [m³] in an allocation area, minimum of
        # available surface water and demand in an allocation area.
        # [m³ s⁻¹] = min([m³] / [s], [m³ s⁻¹])
        sw_abstraction = min(sw_available / dt, sw_demand_vol)

        # fraction of available surface water that can be abstracted at allocation area
        # level
        # [-] = [m³ s⁻¹] / ([m³] / [s])
        frac_abstract_sw = bounded_divide(sw_abstraction, sw_available / dt)
        # fraction of water demand that can be satisfied by available surface water at
        # allocation area level
        # [-] = [m³ s⁻¹] / [m³ s⁻¹]
        frac_allocate_sw = bounded_divide(sw_abstraction, sw_demand_vol)

        # water abstracted from surface water at each river cell (including reservoir
        # locations).
        for j in inds_river[i]
            # [m³ s⁻¹] = [-] * [m³] / [s]
            act_surfacewater_abst_vol[j] +=
                frac_abstract_sw * available_surfacewater[j] / dt
            # [m s⁻¹] = [m³ s⁻¹] / [m²]
            act_surfacewater_abst[j] =
                act_surfacewater_abst_vol[j] / domain.river.parameters.cell_area[j]
        end

        # water allocated to each land cell.
        for j in inds_land[i]
            # [m s⁻¹] = [-] * [m s⁻¹]
            surfacewater_alloc[j] += frac_allocate_sw * surfacewater_demand[j]
        end
    end
end

function available_surface_water!(
    available_surfacewater::Vector{Float64},
    reservoir,
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
            # [m³ s⁻¹]
            external_inflow = reservoir.boundary_conditions.external_inflow[k]
            # [m³] = [m³] * [-]
            available_volume = reservoir.variables.storage[k] * 0.98
            if external_inflow < 0.0
                if available_volume > -external_inflow * dt
                    # [m³] = min([m³ s⁻¹] * [s], [m³])
                    available_volume += external_inflow * dt
                else
                    available_volume = 0.0
                end
            end
            # [m³]
            available_surfacewater[j] = available_volume
            sw_available += available_volume
        else
            # river volume
            # [m³]
            sw_available += available_surfacewater[j]
        end
    end
    return sw_available
end

"Update water allocation for land domain based on local groundwater availability."
function groundwater_allocation_local!(
    model::AllocationLand,
    demand_variables::DemandVariables,
    groundwater_storage::Vector{Float64},
    parameters::LandParameters,
    dt,
)
    (;
        surfacewater_alloc,
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = model.variables
    (; groundwater_demand, total_gross_demand) = demand_variables
    (; area, reservoir_coverage) = parameters

    for i in eachindex(groundwater_demand)
        # groundwater demand based on allocation from surface water.

        # [m s⁻¹]
        groundwater_demand[i] = max(total_gross_demand[i] - surfacewater_alloc[i], 0.0)
        # excluding reservoirs
        if !reservoir_coverage[i]
            # satisfy groundwater demand with available local groundwater volume
            # [m³ s⁻¹] = [m s⁻¹] * [m²]
            groundwater_demand_vol = groundwater_demand[i] * area[i]
            # [m³] = [m³] * [-]
            available_volume = groundwater_storage[i] * 0.75 # limit available groundwater volume
            # [m³] = min([m³ s⁻¹] * [s], [m³])
            abstraction_vol = min(groundwater_demand_vol * dt, available_volume)
            # [m³ s⁻¹] = [m³] / [s]
            actual_groundwater_abstraction_volume = abstraction_vol / dt
            act_groundwater_abst_vol[i] = actual_groundwater_abstraction_volume
            # remaining available groundwater and demand
            # [m³] = max([m³] - [m³], [m³])
            available_groundwater[i] = max(available_volume - abstraction_vol, 0.0)
            # [m s⁻¹] = [m³ s⁻¹] / [m²]
            abstraction = actual_groundwater_abstraction_volume / area[i]
            # [m s⁻¹] = max([m s⁻¹] - [m s⁻¹], [m s⁻¹])
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
function groundwater_allocation_area!(
    model::AllocationLand,
    demand_variables::DemandVariables,
    domain::Domain,
    dt::Float64,
)
    inds_river = domain.river.network.allocation_area_indices
    inds_land = domain.land.network.allocation_area_indices
    (;
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = model.variables

    (; groundwater_demand) = demand_variables
    (; area) = domain.land.parameters

    # loop over allocation areas
    for i in eachindex(inds_river)
        # groundwater demand and availability (allocation area)
        # [m³ s⁻¹]
        gw_demand_vol = 0.0
        # [m³]
        gw_available = 0.0
        for j in inds_land[i]
            # [m³ s⁻¹] = [m s⁻¹] * [m²]
            gw_demand_vol += groundwater_demand[j] * area[j]
            # [m³]
            gw_available += available_groundwater[j]
        end
        # total actual groundwater abstraction [m³] in an allocation area, minimum of
        # available  groundwater and demand in an allocation area.
        # [m³] = min([m³], [m³ s⁻¹] * [s])
        gw_abstraction = min(gw_available, gw_demand_vol * dt)

        # fraction of available groundwater that can be abstracted at allocation area level
        # [-] = [m³] / [m³]
        frac_abstract_gw = bounded_divide(gw_abstraction, gw_available)
        # fraction of water demand that can be satisfied by available groundwater at
        # allocation area level.
        # [-] = [m³] / ([m³ s⁻¹] * [s])
        frac_allocate_gw = bounded_divide(gw_abstraction, gw_demand_vol * dt)

        # water abstracted from groundwater and allocated.
        for j in inds_land[i]
            # [m³ s⁻¹] = [-] * [m³] / [s]
            act_groundwater_abst_vol[j] += frac_abstract_gw * available_groundwater[j] / dt
            # [m s⁻¹] = [m³ s⁻¹] / [m²]
            act_groundwater_abst[j] = act_groundwater_abst_vol[j] / area[j]
            # [m s⁻¹] = [-] * [m s⁻¹]
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
    # [m s⁻¹] = [-] * ([m s⁻¹] + [m s⁻¹])
    @. surfacewater_demand = frac_sw_used * (nonirri_demand_gross + irri_demand_gross)

    # local surface water demand and allocation (river, excluding reservoirs)
    surface_water_allocation_local!(model, demand.variables, river, domain.land, dt)
    # surface water demand and allocation for areas
    surface_water_allocation_area!(model, demand.variables, river, domain, dt)

    # [m³ s⁻¹]
    @. abstraction = act_surfacewater_abst_vol

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
        demand.variables,
        groundwater_storage(routing.subsurface_flow),
        domain.land.parameters,
        dt,
    )
    # groundwater demand and allocation for areas
    groundwater_allocation_area!(model, demand.variables, domain)

    # irrigation allocation
    for i in eachindex(total_alloc)
        # [m s⁻¹] = [m s⁻¹] + [m s⁻¹]
        total_alloc[i] = groundwater_alloc[i] + surfacewater_alloc[i]
        # [-] = [m s⁻¹] / [m s⁻¹]
        frac_irri = bounded_divide(irri_demand_gross[i], total_gross_demand[i])
        # [m s⁻¹] = [-] * [m s⁻¹]
        irri_alloc[i] = frac_irri * total_alloc[i]
        # [m s⁻¹] = [m s⁻¹] - [m s⁻¹]
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
            # [m s⁻¹] = [m s⁻¹]
            river.allocation.variables.nonirri_returnflow[k] = nonirri_returnflow[i]
            nonirri_returnflow[i] = 0.0
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
function update_water_demand!(model::Demand, soil::SbmSoilModel, dt::Number)
    (; nonpaddy, paddy, domestic, industry, livestock) = model

    return_flow_fraction!(industry)
    return_flow_fraction!(domestic)
    return_flow_fraction!(livestock)

    update_demand_gross!(nonpaddy, soil, dt)
    update_demand_gross!(paddy, dt)
    update_demand_gross!(model)

    return nothing
end

update_water_demand!(model::NoDemand, soil::SbmSoilModel, dt::Number) = nothing
