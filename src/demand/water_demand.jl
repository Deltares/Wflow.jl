abstract type AbstractIrrigationModel{T} end
abstract type AbstractAllocationModel{T} end
abstract type AbstractDemandModel end

struct NoIrrigationPaddy{T} <: AbstractIrrigationModel{T} end
struct NoIrrigationNonPaddy{T} <: AbstractIrrigationModel{T} end
struct NoNonIrrigationDemand <: AbstractDemandModel end
struct NoAllocationLand{T} <: AbstractAllocationModel{T} end
struct NoAllocationRiver{T} <: AbstractAllocationModel{T} end

"Struct to store non-irrigation water demand variables"
@get_units @grid_loc @with_kw struct NonIrrigationDemandVariables{T}
    returnflow::Vector{T}                   # return flow [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
end

"Struct to store prescribed water demand variables"
@get_units @grid_loc @with_kw struct PrescibedDemand{T}
    demand_gross::Vector{T}     # gross water demand [mm Δt⁻¹]
    demand_net::Vector{T}       # net water demand [mm Δt⁻¹]
end

"Non-irrigation water demand model"
@with_kw struct NonIrrigationDemand{T} <: AbstractDemandModel
    demand::PrescibedDemand{T}
    variables::NonIrrigationDemandVariables{T}
end

# wrapper methods
get_demand_gross(model::NonIrrigationDemand) = model.demand.demand_gross
get_demand_gross(model::NoNonIrrigationDemand) = 0.0

"Initialize non-irrigation water demand model for a water use `sector`"
function NonIrrigationDemand(dataset, config, indices, dt, sector)
    lens = lens_input_parameter("land~$(sector)__gross_water_demand_volume_flux")
    demand_gross =
        ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float) .*
        (dt / basetimestep)
    lens = lens_input_parameter("land~$(sector)__net_water_demand_volume_flux")
    demand_net =
        ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float) .*
        (dt / basetimestep)
    n = length(indices)
    returnflow_f = return_flow_fraction.(demand_gross, demand_net)

    demand = PrescibedDemand{Float}(; demand_gross, demand_net)
    vars = NonIrrigationDemandVariables{Float}(;
        returnflow_fraction = returnflow_f,
        returnflow = fill(Float(0), n),
    )
    non_irrigation_demand = NonIrrigationDemand{Float}(; demand, variables = vars)

    return non_irrigation_demand
end

"Struct to store non-paddy irrigation model variables"
@get_units @grid_loc @with_kw struct NonPaddyVariables{T}
    demand_gross::Vector{T}     # irrigation gross demand [mm Δt⁻¹] 
end

"Struct to store non-paddy irrigation model parameters"
@get_units @grid_loc @with_kw struct NonPaddyParameters{T}
    irrigation_efficiency::Vector{T} | "-"      # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{T}          # maximum irrigation rate [mm Δt⁻¹]
    irrigation_areas::Vector{Bool} | "-"        # irrigation areas [-]
    irrigation_trigger::Vector{Bool} | "-"      # irrigation on or off [-]end
end

"Non-paddy (other crops than flooded rice) irrigation model"
@with_kw struct NonPaddy{T} <: AbstractIrrigationModel{T}
    parameters::NonPaddyParameters{T}
    variables::NonPaddyVariables{T}
end

"Initialize non-paddy irrigation model"
function NonPaddy(dataset, config, indices, dt)
    lens = lens_input_parameter("land~irrigated-non-paddy__irrigation_efficiency")
    efficiency = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float)

    lens = lens_input_parameter("land~irrigated-non-paddy_area__number")
    areas = ncread(
        dataset,
        config,
        lens;
        sel = indices,
        defaults = 1,
        optional = false,
        type = Int,
    )
    lens = lens_input_parameter("land~irrigated-non-paddy__irrigation_trigger_flag")
    irrigation_trigger = ncread(
        dataset,
        config,
        lens;
        sel = indices,
        defaults = 1,
        optional = false,
        type = Bool,
    )
    lens = lens_input_parameter("land~irrigated-non-paddy__max_irrigation_rate")
    max_irri_rate =
        ncread(dataset, config, lens; sel = indices, defaults = 25.0, type = Float) .*
        (dt / basetimestep)

    params = NonPaddyParameters{Float}(;
        maximum_irrigation_rate = max_irri_rate,
        irrigation_efficiency = efficiency,
        irrigation_areas = areas,
        irrigation_trigger,
    )
    vars = NonPaddyVariables{Float}(; demand_gross = fill(mv, length(indices)))

    nonpaddy = NonPaddy{Float}(; variables = vars, parameters = params)

    return nonpaddy
end

# wrapper methods
get_demand_gross(model::NonPaddy) = model.variables.demand_gross
get_demand_gross(model::NoIrrigationNonPaddy) = 0.0

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
    (; h3, n_unsatlayers, zi, ustorelayerdepth, f_infiltration_reduction) = soil.variables
    (;
        irrigation_areas,
        irrigation_trigger,
        maximum_irrigation_rate,
        irrigation_efficiency,
    ) = model.parameters
    rootingdepth = get_rootingdepth(soil)

    for i in eachindex(irrigation_areas)
        if irrigation_areas[i] && irrigation_trigger[i]
            usl = set_layerthickness(zi[i], sumlayers[i], act_thickl[i])
            irri_dem_gross = 0.0
            for k in 1:n_unsatlayers[i]
                # compute water demand only for root zone through root fraction per layer
                rootfrac = min(1.0, (max(0.0, rootingdepth[i] - sumlayers[i][k]) / usl[k]))
                # vwc_f and vwc_h3 can be precalculated.
                vwc_fc = vwc_brooks_corey(-100.0, hb[i], theta_s[i], theta_r[i], c[i][k])
                vwc_h3 = vwc_brooks_corey(h3[i], hb[i], theta_s[i], theta_r[i], c[i][k])
                depletion =
                    (vwc_fc * usl[k]) - (ustorelayerdepth[i][k] + theta_r[i] * usl[k])
                depletion *= rootfrac
                raw = (vwc_fc - vwc_h3) * usl[k] # readily available water
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
@get_units @grid_loc @with_kw struct PaddyVariables{T}
    demand_gross::Vector{T}     # irrigation gross demand [mm Δt⁻¹]
    h::Vector{T} | "mm"         # actual water depth in rice field [mm]
    evaporation::Vector{T}      # evaporation rate [mm Δt⁻¹] 
end

"Struct to store paddy irrigation model parameters"
@get_units @grid_loc @with_kw struct PaddyParameters{T}
    irrigation_efficiency::Vector{T} | "-"      # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{T}          # maximum irrigation rate [mm Δt⁻¹]
    irrigation_areas::Vector{Bool} | "-"        # irrigation areas [-]
    irrigation_trigger::Vector{Bool} | "-"      # irrigation on or off [-]
    h_min::Vector{T} | "mm"                     # minimum required water depth in the irrigated rice field [mm]
    h_opt::Vector{T} | "mm"                     # optimal water depth in the irrigated rice fields [mm]
    h_max::Vector{T} | "mm"                     # water depth when rice field starts spilling water (overflow) [mm]
end

"Paddy (flooded rice) irrigation model"
@with_kw struct Paddy{T} <: AbstractIrrigationModel{T}
    parameters::PaddyParameters{T}
    variables::PaddyVariables{T}
end

"Initialize paddy irrigation model"
function Paddy(dataset, config, indices, dt)
    lens = lens_input_parameter("land~irrigated-paddy__min_depth")
    h_min = ncread(dataset, config, lens; sel = indices, defaults = 20.0, type = Float)

    lens = lens_input_parameter("land~irrigated-paddy__optimal_depth")
    h_opt = ncread(dataset, config, lens; sel = indices, defaults = 50.0, type = Float)

    lens = lens_input_parameter("land~irrigated-paddy__max_depth")
    h_max = ncread(dataset, config, lens; sel = indices, defaults = 80.0, type = Float)

    lens = lens_input_parameter("land~irrigated-paddy__irrigation_efficiency")
    efficiency = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float)

    lens = lens_input_parameter("land~irrigated-paddy_area__number")
    areas = ncread(dataset, config, lens; sel = indices, optional = false, type = Bool)

    lens = lens_input_parameter("land~irrigated-paddy__irrigation_trigger_flag")
    irrigation_trigger =
        ncread(dataset, config, lens; sel = indices, optional = false, type = Bool)
    lens = lens_input_parameter("land~irrigate-paddy__max_irrigation_rate")
    max_irri_rate =
        ncread(dataset, config, lens; sel = indices, defaults = 25.0, type = Float) .*
        (dt / basetimestep)
    n = length(indices)
    params = PaddyParameters{Float}(;
        irrigation_efficiency = efficiency,
        maximum_irrigation_rate = max_irri_rate,
        irrigation_trigger,
        h_min,
        h_max,
        h_opt,
        irrigation_areas = areas,
    )
    vars = PaddyVariables{Float}(;
        demand_gross = fill(mv, n),
        h = fill(0.0, n),
        evaporation = fill(0.0, n),
    )
    paddy = Paddy{Float}(; parameters = params, variables = vars)
    return paddy
end

# wrapper methods
get_water_depth(model::Paddy) = model.variables.h
get_water_depth(model::NoIrrigationPaddy) = 0.0
get_demand_gross(model::Paddy) = model.variables.demand_gross
get_demand_gross(model::NoIrrigationPaddy) = 0.0

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
get_evaporation(model::NoIrrigationPaddy) = 0.0
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
@get_units @grid_loc @with_kw struct DemandVariables{T}
    irri_demand_gross::Vector{T}        # irrigation gross demand [mm Δt⁻¹]
    nonirri_demand_gross::Vector{T}     # non-irrigation gross demand [mm Δt⁻¹]
    total_gross_demand::Vector{T}       # total gross demand [mm Δt⁻¹]
    surfacewater_demand::Vector{T}      # demand from surface water [mm Δt⁻¹]
    groundwater_demand::Vector{T}       # demand from groundwater [mm Δt⁻¹]
end

"Initialize water demand variables"
function DemandVariables(T::Type{<:AbstractFloat}, n::Int)
    return DemandVariables{T}(;
        irri_demand_gross = zeros(T, n),
        nonirri_demand_gross = zeros(T, n),
        total_gross_demand = zeros(T, n),
        surfacewater_demand = zeros(T, n),
        groundwater_demand = zeros(T, n),
    )
end

"Water demand model"
@with_kw struct Demand{D, I, L, P, NP, V} <: AbstractDemandModel
    domestic::D
    industry::I
    livestock::L
    paddy::P
    nonpaddy::NP
    variables::V
end

@with_kw struct NoDemand{T} <: AbstractDemandModel
    domestic::NoNonIrrigationDemand = NoNonIrrigationDemand()
    industry::NoNonIrrigationDemand = NoNonIrrigationDemand()
    livestock::NoNonIrrigationDemand = NoNonIrrigationDemand()
    paddy::NoIrrigationPaddy{T} = NoIrrigationPaddy{T}()
    nonpaddy::NoIrrigationNonPaddy{T} = NoIrrigationNonPaddy{T}()
end

"Initialize water demand model"
function Demand(dataset, config, indices, dt)
    domestic = if get(config.model.water_demand, "domestic", false)
        NonIrrigationDemand(dataset, config, indices, dt, "domestic")
    else
        NoNonIrrigationDemand()
    end
    industry = if get(config.model.water_demand, "industry", false)
        NonIrrigationDemand(dataset, config, indices, dt, "industry")
    else
        NoNonIrrigationDemand()
    end
    livestock = if get(config.model.water_demand, "livestock", false)
        NonIrrigationDemand(dataset, config, indices, dt, "livestock")
    else
        NoNonIrrigationDemand()
    end
    paddy = if get(config.model.water_demand, "paddy", false)
        Paddy(dataset, config, indices, dt)
    else
        NoIrrigationPaddy{Float}()
    end
    nonpaddy = if get(config.model.water_demand, "nonpaddy", false)
        NonPaddy(dataset, config, indices, dt)
    else
        NoIrrigationNonPaddy{Float}()
    end

    n = length(indices)
    vars = DemandVariables(Float, n)
    demand = Demand(; domestic, industry, livestock, paddy, nonpaddy, variables = vars)
    return demand
end

"Struct to store river allocation model variables"
@get_units @grid_loc @with_kw struct AllocationRiverVariables{T}
    act_surfacewater_abst::Vector{T}                        # actual surface water abstraction [mm Δt⁻¹]
    act_surfacewater_abst_vol::Vector{T} | "m3 dt-1"        # actual surface water abstraction [m³ Δt⁻¹]
    available_surfacewater::Vector{T} | "m3"                # available surface water [m³]
    nonirri_returnflow::Vector{T}                           # return flow from non irrigation [mm Δt⁻¹] 
end

"Initialize river allocation model variables"
function AllocationRiverVariables(T::Type{<:AbstractFloat}, n::Int)
    return AllocationRiverVariables{T}(;
        act_surfacewater_abst = zeros(T, n),
        act_surfacewater_abst_vol = zeros(T, n),
        available_surfacewater = zeros(T, n),
        nonirri_returnflow = zeros(T, n),
    )
end

"River allocation model"
@with_kw struct AllocationRiver{T} <: AbstractAllocationModel{T}
    variables::AllocationRiverVariables{T}
end

get_nonirrigation_returnflow(model::AllocationRiver) = model.variables.nonirri_returnflow
get_nonirrigation_returnflow(model::NoAllocationRiver) = 0.0

"Initialize water allocation for the river domain"
function AllocationRiver(n)
    vars = AllocationRiverVariables(Float, n)
    allocation = AllocationRiver{Float}(; variables = vars)
    return allocation
end

"Struct to store land allocation allocation model parameters"
@get_units @grid_loc @with_kw struct AllocationLandParameters{T}
    frac_sw_used::Vector{T} | "-"       # fraction surface water used [-]
    areas::Vector{Int} | "-"            # allocation areas [-]
end

"Struct to store land allocation model variables"
@get_units @grid_loc @with_kw struct AllocationLandVariables{T}
    surfacewater_alloc::Vector{T}                       # allocation from surface water [mm Δt⁻¹]
    act_groundwater_abst::Vector{T}                     # actual groundwater abstraction [mm Δt⁻¹]
    act_groundwater_abst_vol::Vector{T} | "m3 dt-1"     # actual groundwater abstraction [m³ Δt⁻¹]
    available_groundwater::Vector{T} | "m3"             # available groundwater [m³]
    groundwater_alloc::Vector{T}                        # allocation from groundwater [mm Δt⁻¹]
    irri_alloc::Vector{T}                               # allocated water for irrigation [mm Δt⁻¹]
    nonirri_alloc::Vector{T}                            # allocated water for non-irrigation [mm Δt⁻¹]
    total_alloc::Vector{T}                              # total allocated water [mm Δt⁻¹]
    nonirri_returnflow::Vector{T}                       # return flow from non irrigation [mm Δt⁻¹]
end

"Initialize land allocation model variables"
function AllocationLandVariables(T::Type{<:AbstractFloat}, n::Int)
    return AllocationLandVariables{T}(;
        surfacewater_alloc = zeros(T, n),
        act_groundwater_abst = zeros(T, n),
        act_groundwater_abst_vol = zeros(T, n),
        available_groundwater = zeros(T, n),
        groundwater_alloc = zeros(T, n),
        irri_alloc = zeros(T, n),
        nonirri_alloc = zeros(T, n),
        total_alloc = zeros(T, n),
        nonirri_returnflow = zeros(T, n),
    )
end

"Land allocation model"
@with_kw struct AllocationLand{T} <: AbstractAllocationModel{T}
    parameters::AllocationLandParameters{T}
    variables::AllocationLandVariables{T}
end

"Initialize water allocation for the land domain"
function AllocationLand(dataset, config, indices)
    lens = lens_input_parameter("land_surface_water__withdrawal_fraction")
    frac_sw_used = ncread(dataset, config, lens; sel = indices, defaults = 1, type = Float)

    lens = lens_input_parameter("land_water_allocation_area__number")
    areas = ncread(dataset, config, lens; sel = indices, defaults = 1, type = Int)

    n = length(indices)

    params = AllocationLandParameters(; areas = areas, frac_sw_used = frac_sw_used)
    vars = AllocationLandVariables(Float, n)
    allocation = AllocationLand(; parameters = params, variables = vars)
    return allocation
end

# wrapper methods
get_irrigation_allocated(model::AllocationLand) = model.variables.irri_alloc
get_irrigation_allocated(model::NoAllocationLand) = 0.0
get_nonirrigation_returnflow(model::AllocationLand) = model.variables.nonirri_returnflow
get_nonirrigation_returnflow(model::NoAllocationLand) = 0.0

"""
Return return flow fraction based on gross water demand `demand_gross` and net water demand
`demand_net`
"""
function return_flow_fraction(demand_gross, demand_net)
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
function surface_water_allocation_local!(model::AllocationLand, demand, river, network, dt)
    (; surfacewater_alloc) = model.variables
    (; surfacewater_demand) = demand.variables
    (; act_surfacewater_abst_vol, act_surfacewater_abst, available_surfacewater) =
        river.allocation.variables
    (; inflow) = river.boundary_conditions
    (; volume) = river.variables
    # maps from the land domain to the internal river domain (linear index), excluding water bodies
    index_river = network.land.river_inds_excl_waterbody
    for i in eachindex(surfacewater_demand)
        if index_river[i] > 0.0
            # the available volume is limited by a fixed scaling factor of 0.8 to prevent
            # rivers completely drying out. check for abstraction through inflow (external
            # negative inflow) and adjust available volume.
            if inflow[index_river[i]] < 0.0
                river_inflow = inflow[index_river[i]] * dt
                available_volume = max(volume[index_river[i]] * 0.80 + river_inflow, 0.0)
            else
                available_volume = volume[index_river[i]] * 0.80
            end
            # satisfy surface water demand with available local river volume
            surfacewater_demand_vol = surfacewater_demand[i] * 0.001 * network.land.area[i]
            abstraction_vol = min(surfacewater_demand_vol, available_volume)
            act_surfacewater_abst_vol[index_river[i]] = abstraction_vol
            # remaining available surface water and demand 
            available_surfacewater[index_river[i]] =
                max(available_volume - abstraction_vol, 0.0)
            abstraction = (abstraction_vol / network.land.area[i]) * 1000.0
            surfacewater_demand[i] = max(surfacewater_demand[i] - abstraction, 0.0)
            # update actual abstraction from river and surface water allocation (land cell)
            act_surfacewater_abst[index_river[i]] = abstraction
            surfacewater_alloc[i] = abstraction
        end
    end
    return nothing
end

"""
Update water allocation for river and land domains based on surface water (river)
availability for allocation areas.
"""
function surface_water_allocation_area!(model::AllocationLand, demand, river, network)
    inds_river = network.river.allocation_area_indices
    inds_land = network.land.allocation_area_indices
    inds_reservoir = network.river.reservoir_indices
    inds_lake = network.river.lake_indices

    (; available_surfacewater, act_surfacewater_abst_vol, act_surfacewater_abst) =
        river.allocation.variables
    (; surfacewater_alloc) = model.variables
    (; surfacewater_demand) = demand.variables
    (; reservoir, lake) = river.boundary_conditions

    # loop over allocation areas
    for i in eachindex(inds_river)
        # surface water demand (allocation area)
        sw_demand_vol = 0.0
        for j in inds_land[i]
            sw_demand_vol += surfacewater_demand[j] * 0.001 * network.land.area[j]
        end
        # surface water availability (allocation area)
        sw_available = 0.0
        for j in inds_river[i]
            if inds_reservoir[j] > 0
                # for reservoir locations use reservoir volume
                k = inds_reservoir[j]
                available_surfacewater[j] = reservoir.volume[k] * 0.98 # limit available reservoir volume
                sw_available += available_surfacewater[j]
            elseif inds_lake[j] > 0
                # for lake locations use lake volume
                k = inds_lake[j]
                available_surfacewater[j] = lake.storage[k] * 0.98 # limit available lake volume
                sw_available += available_surfacewater[j]

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

        # water abstracted from surface water at each river cell (including reservoir and
        # lake locations).
        for j in inds_river[i]
            act_surfacewater_abst_vol[j] += frac_abstract_sw * available_surfacewater[j]
            act_surfacewater_abst[j] =
                (act_surfacewater_abst_vol[j] / network.river.cell_area[j]) * 1000.0
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
    demand,
    groundwater_volume,
    network,
)
    (;
        surfacewater_alloc,
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = model.variables
    (; groundwater_demand, total_gross_demand) = demand.variables

    for i in eachindex(groundwater_demand)
        # groundwater demand based on allocation from surface water.
        groundwater_demand[i] = max(total_gross_demand[i] - surfacewater_alloc[i], 0.0)
        # excluding water bodies
        if !network.waterbody[i]
            # satisfy groundwater demand with available local groundwater volume
            groundwater_demand_vol = groundwater_demand[i] * 0.001 * network.area[i]
            available_volume = groundwater_volume[i] * 0.75 # limit available groundwater volume
            abstraction_vol = min(groundwater_demand_vol, available_volume)
            act_groundwater_abst_vol[i] = abstraction_vol
            # remaining available groundwater and demand 
            available_groundwater[i] = max(available_volume - abstraction_vol, 0.0)
            abstraction = (abstraction_vol / network.area[i]) * 1000.0
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
function groundwater_allocation_area!(model::AllocationLand, demand, network)
    inds_river = network.river.allocation_area_indices
    inds_land = network.land.allocation_area_indices
    (;
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = model.variables

    (; groundwater_demand) = demand.variables

    # loop over allocation areas
    for i in eachindex(inds_river)
        # groundwater demand and availability (allocation area)
        gw_demand_vol = 0.0
        gw_available = 0.0
        for j in inds_land[i]
            gw_demand_vol += groundwater_demand[j] * 0.001 * network.land.area[j]
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
            act_groundwater_abst[j] =
                1000.0 * (act_groundwater_abst_vol[j] / network.land.area[j])
            groundwater_alloc[j] += frac_allocate_gw * groundwater_demand[j]
        end
    end
    return nothing
end

"Return and update non-irrigation sector (domestic, livestock, industry) return flow"
function return_flow(model::NonIrrigationDemand, nonirri_demand_gross, nonirri_alloc)
    for i in eachindex(model.variables.returnflow)
        frac = bounded_divide(model.demand.demand_gross[i], nonirri_demand_gross[i])
        allocate = frac * nonirri_alloc[i]
        model.variables.returnflow[i] = model.variables.returnflow_fraction[i] * allocate
    end
    return model.variables.returnflow
end

# return zero (return flow) if non-irrigation sector is not defined
return_flow(model::NoNonIrrigationDemand, nonirri_demand_gross, nonirri_alloc) = 0.0

# wrapper methods
groundwater_volume(model::LateralSSF) = model.variables.volume
groundwater_volume(model) = model.flow.aquifer.variables.volume

"""
    update_water_allocation!((model::AllocationLand, demand, lateral, network, dt)

Update water allocation for the land domain `AllocationLand` and water allocation for the
river domain (part of `lateral`) based on the water `demand` model for a single timestep.
First, surface water abstraction is computed to satisfy local water demand (non-irrigation
and irrigation), and then updated (including lakes and reservoirs) to satisfy the remaining
water demand for allocation areas. Then groundwater abstraction is computed to satisfy the
remaining local water demand, and then updated to satisfy the remaining water demand for
allocation areas. Finally, non-irrigation return flows are updated.
"""
function update_water_allocation!(model::AllocationLand, demand, lateral, network, dt)
    river = lateral.river
    index_river = network.land.river_inds_excl_waterbody
    inds_reservoir = network.reservoir.river_indices
    inds_lake = network.lake.river_indices
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
    (; abstraction, reservoir, lake) = river.boundary_conditions

    surfacewater_alloc .= 0.0
    act_surfacewater_abst .= 0.0
    act_surfacewater_abst_vol .= 0.0
    # total surface water demand for each land cell
    @. surfacewater_demand =
        frac_sw_used * nonirri_demand_gross + frac_sw_used * irri_demand_gross

    # local surface water demand and allocation (river, excluding reservoirs and lakes)
    surface_water_allocation_local!(model, demand, river, network, dt)
    # surface water demand and allocation for areas
    surface_water_allocation_area!(model, demand, river, network)

    @. abstraction = act_surfacewater_abst_vol / dt

    # for reservoir and lake locations set river abstraction at zero and abstract volume
    # from reservoir and lake, including an update of lake waterlevel
    if !isnothing(reservoir)
        @. abstraction[inds_reservoir] = 0.0
        @. reservoir.volume -= act_surfacewater_abst_vol[inds_reservoir]
    elseif !isnothing(lake)
        @. abstraction[inds_lake] = 0.0
        @. lake.storage -= act_surfacewater_abst_vol[inds_lake]
        @. lake.waterlevel = waterlevel(lake.storfunc, lake.area, lake.storage, lake.sh)
    end

    groundwater_alloc .= 0.0
    act_groundwater_abst_vol .= 0.0
    act_groundwater_abst .= 0.0
    # local groundwater demand and allocation
    groundwater_allocation_local!(
        model,
        demand,
        groundwater_volume(lateral.subsurface),
        network.land,
    )
    # groundwater demand and allocation for areas
    groundwater_allocation_area!(model, demand, network)

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
        if index_river[i] > 0.0
            k = index_river[i]
            river.allocation.variables.nonirri_returnflow[k] = nonirri_returnflow[i]
            nonirri_returnflow[i] = 0.0
        else
            nonirri_returnflow[i] = nonirri_returnflow[i]
        end
    end
end
update_water_allocation!(model::NoAllocationLand, demand, lateral, network, dt) = nothing

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
    update_water_demand!(model::Demand, soil)

Update the return flow fraction `returnflow_fraction` of `industry`, `domestic` and
`livestock`, gross water demand `demand_gross` of `paddy` and `nonpaddy` models, and the
total gross water demand, total irrigation gross water demand and total non-irrigation gross
water demand as part of the water `demand` model.
"""
function update_water_demand!(model::Demand, soil)
    (; nonpaddy, paddy, domestic, industry, livestock) = model

    return_flow_fraction!(industry)
    return_flow_fraction!(domestic)
    return_flow_fraction!(livestock)

    update_demand_gross!(nonpaddy, soil)
    update_demand_gross!(paddy)
    update_demand_gross!(model)

    return nothing
end
update_water_demand!(model::NoDemand, soil) = nothing