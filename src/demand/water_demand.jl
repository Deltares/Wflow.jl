abstract type AbstractIrrigationModel{T} end
abstract type AbstractAllocationModel{T} end
abstract type AbstractDemandModel end

struct NoIrrigationPaddy{T} <: AbstractIrrigationModel{T} end
struct NoIrrigationNonPaddy{T} <: AbstractIrrigationModel{T} end
struct NoNonIrrigationDemand <: AbstractDemandModel end
struct NoAllocationLand{T} <: AbstractAllocationModel{T} end

@get_units @grid_loc @with_kw struct NonIrrigationDemandVariables{T}
    returnflow::Vector{T}                   # return flow [mm Δt⁻¹]
    returnflow_fraction::Vector{T} | "-"    # return flow fraction [-]
end

@get_units @grid_loc @with_kw struct PrescibedDemand{T}
    demand_gross::Vector{T}                 # gross water demand [mm Δt⁻¹]
    demand_net::Vector{T}                   # net water demand [mm Δt⁻¹]
end

@get_units @grid_loc @with_kw struct NonIrrigationDemand{T} <: AbstractDemandModel
    demand::PrescibedDemand{T} | "-"
    variables::NonIrrigationDemandVariables{T} | "-"
end

get_demand_gross(model::NonIrrigationDemand) = model.demand.demand_gross
get_demand_gross(model::NoNonIrrigationDemand) = 0.0

@get_units @grid_loc @with_kw struct NonPaddyVariables{T}
    demand_gross::Vector{T}                     # irrigation gross demand [mm Δt⁻¹] 
end

@get_units @grid_loc @with_kw struct NonPaddyParameters{T}
    irrigation_efficiency::Vector{T} | "-"      # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{T}          # maximum irrigation rate [mm Δt⁻¹]
    irrigation_areas::Vector{Bool} | "-"        # irrigation areas [-]
    irrigation_trigger::Vector{Bool} | "-"      # irrigation on or off [-]end
end

@get_units @grid_loc @with_kw struct NonPaddy{T} <: AbstractIrrigationModel{T}
    parameters::NonPaddyParameters{T} | "-"
    variables::NonPaddyVariables{T} | "-"
end

get_demand_gross(model::NonPaddy) = model.variables.demand_gross
get_demand_gross(model::NoIrrigationNonPaddy) = 0.0

function update_demand_gross!(nonpaddy::NonPaddy, soil::SbmSoilModel)
    (; hb, theta_s, theta_r, c, sumlayers, act_thickl, pathfrac, infiltcapsoil) =
        soil.parameters
    (; h3, n_unsatlayers, zi, ustorelayerdepth, f_infiltration_reduction) = soil.variables
    (;
        irrigation_areas,
        irrigation_trigger,
        maximum_irrigation_rate,
        irrigation_efficiency,
    ) = nonpaddy.parameters
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
                    nonpaddy.variables.demand_gross[i] == maximum_irrigation_rate[i]
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
        nonpaddy.variables.demand_gross[i] = irri_dem_gross
    end
end
update_demand_gross!(nonpaddy::NoIrrigationNonPaddy, soil::SbmSoilModel) = nothing

@get_units @grid_loc @with_kw struct PaddyVariables{T}
    demand_gross::Vector{T}                     # irrigation gross demand [mm Δt⁻¹]
    h::Vector{T} | "mm"                         # actual water depth in rice field [mm]
    evaporation::Vector{T}                      # evaporation rate [mm Δt⁻¹] 
end

@get_units @grid_loc @with_kw struct PaddyParameters{T}
    irrigation_efficiency::Vector{T} | "-"      # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{T}          # maximum irrigation rate [mm Δt⁻¹]
    irrigation_areas::Vector{Bool} | "-"        # irrigation areas [-]
    irrigation_trigger::Vector{Bool} | "-"      # irrigation on or off [-]
    h_min::Vector{T} | "mm"                     # minimum required water depth in the irrigated rice field [mm]
    h_opt::Vector{T} | "mm"                     # optimal water depth in the irrigated rice fields [mm]
    h_max::Vector{T} | "mm"                     # water depth when rice field starts spilling water (overflow) [mm]
end

@get_units @grid_loc @with_kw struct Paddy{T} <: AbstractIrrigationModel{T}
    parameters::PaddyParameters{T} | "-"
    variables::PaddyVariables{T} | "-"
end

get_water_depth(model::Paddy) = model.variables.h
get_water_depth(model::NoIrrigationPaddy) = 0.0
get_demand_gross(model::Paddy) = model.variables.demand_gross
get_demand_gross(model::NoIrrigationPaddy) = 0.0

function evaporation!(model::Paddy, potential_evaporation)
    for i in eachindex(potential_evaporation)
        if model.parameters.irrigation_areas[i]
            evaporation = min(model.variables.h[i], potential_evaporation[i])
            model.variables.h[i] -= evaporation
            model.variables.evaporation[i] = evaporation
        end
    end
end
evaporation!(model::NoIrrigationPaddy, potential_evaporation) = nothing

get_evaporation(model::NoIrrigationPaddy) = 0.0
get_evaporation(model::Paddy) = model.variables.evaporation

function update_runoff!(model::Paddy, runoff)
    for i in eachindex(model.parameters.irrigation_areas)
        if model.parameters.irrigation_areas[i]
            paddy_runoff = max(runoff[i] - model.parameters.h_max[i], 0.0)
            model.variables.h[i] = runoff[i] - paddy_runoff
            runoff[i] = paddy_runoff
        end
    end
    return runoff
end

function update_runoff!(model::NoIrrigationPaddy, runoff)
    return runoff
end

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
            if h[i] < h_min[i]
                irr_depth_paddy = h_opt[i] - h[i]
            elseif h[i] < h_opt[i] && max_irri_rate_applied # continue irrigation
                irr_depth_paddy = h_opt[i] - h[i]
            else
                irr_depth_paddy = 0.0
            end
            irri_dem_gross = irr_depth_paddy / irrigation_efficiency[i]
            # limit irrigation demand to the maximum irrigation rate
            irri_dem_gross = min(irri_dem_gross, maximum_irrigation_rate[i])
        else
            irri_dem_gross = 0.0
        end
        demand_gross[i] = irri_dem_gross
    end
end

update_demand_gross!(paddy::NoIrrigationPaddy) = nothing

@get_units @grid_loc @with_kw struct DemandVariables{T}
    irri_demand_gross::Vector{T}                        # irrigation gross demand [mm Δt⁻¹]
    nonirri_demand_gross::Vector{T}                     # non-irrigation gross demand [mm Δt⁻¹]
    total_gross_demand::Vector{T}                       # total gross demand [mm Δt⁻¹]
    surfacewater_demand::Vector{T}                      # demand from surface water [mm Δt⁻¹]
    groundwater_demand::Vector{T}                       # demand from groundwater [mm Δt⁻¹]
end

@get_units @grid_loc @with_kw struct Demand{D, I, L, P, NP, V} <: AbstractDemandModel
    domestic::D | "-" | "none"
    industry::I | "-" | "none"
    livestock::L | "-" | "none"
    paddy::P | "-" | "none"
    nonpaddy::NP | "-" | "none"
    variables::V | "-" | "none"
end

@get_units @grid_loc @with_kw struct NoDemand{T} <: AbstractDemandModel
    domestic::NoNonIrrigationDemand = NoNonIrrigationDemand() | "-" | "none"
    industry::NoNonIrrigationDemand = NoNonIrrigationDemand() | "-" | "none"
    livestock::NoNonIrrigationDemand = NoNonIrrigationDemand() | "-" | "none"
    paddy::NoIrrigationPaddy{T} = NoIrrigationPaddy{T}() | "-" | "none"
    nonpaddy::NoIrrigationNonPaddy{T} = NoIrrigationNonPaddy{T}() | "-" | "none"
end

@get_units @grid_loc @with_kw struct AllocationRiverVariables{T}
    act_surfacewater_abst::Vector{T}                    # actual surface water abstraction [mm Δt⁻¹]
    act_surfacewater_abst_vol::Vector{T} | "m3 dt-1"    # actual surface water abstraction [m³ Δt⁻¹]
    available_surfacewater::Vector{T} | "m3"            # available surface water [m³]
    nonirri_returnflow::Vector{T}                       # return flow from non irrigation [mm Δt⁻¹] 
end

@get_units @grid_loc @with_kw struct AllocationRiver{T} <: AbstractAllocationModel{T}
    variables::AllocationRiverVariables{T}
end

@get_units @grid_loc @with_kw struct AllocationLandParameters{T}
    frac_sw_used::Vector{T} | "-"                       # fraction surface water used [-]
    areas::Vector{Int} | "-"                            # allocation areas [-]
end

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

@get_units @grid_loc @with_kw struct AllocationLand{T} <: AbstractAllocationModel{T}
    parameters::AllocationLandParameters{T}
    variables::AllocationLandVariables{T}
end

get_irrigation_allocated(model::AllocationLand) = model.variables.irri_alloc
get_irrigation_allocated(model::NoAllocationLand) = 0.0

"Return return flow fraction based on gross water demand `demand_gross` and net water demand `demand_net`"
function return_flow_fraction(demand_gross, demand_net)
    fraction = bounded_divide(demand_net, demand_gross)
    returnflow_fraction = 1.0 - fraction
    return returnflow_fraction
end

function Demand(nc, config, inds, dt)
    if get(config.model.water_demand, "domestic", false)
        domestic = NonIrrigationDemand(nc, config, inds, dt; sector = "domestic")
    else
        domestic = NoNonIrrigationDemand()
    end
    if get(config.model.water_demand, "industry", false)
        industry = NonIrrigationDemand(nc, config, inds, dt; sector = "industry")
    else
        industry = NoNonIrrigationDemand()
    end
    if get(config.model.water_demand, "livestock", false)
        livestock = NonIrrigationDemand(nc, config, inds, dt; sector = "livestock")
    else
        livestock = NoNonIrrigationDemand()
    end
    if get(config.model.water_demand, "paddy", false)
        paddy = Paddy(nc, config, inds, dt)
    else
        paddy = NoIrrigationPaddy{Float}()
    end
    if get(config.model.water_demand, "nonpaddy", false)
        nonpaddy = NonPaddy(nc, config, inds, dt)
    else
        nonpaddy = NoIrrigationNonPaddy{Float}()
    end

    n = length(inds)
    vars = DemandVariables{Float}(;
        irri_demand_gross = zeros(Float, n),
        nonirri_demand_gross = zeros(Float, n),
        total_gross_demand = zeros(Float, n),
        surfacewater_demand = zeros(Float, n),
        groundwater_demand = zeros(Float, n),
    )

    demand = Demand(;
        domestic = domestic,
        industry = industry,
        livestock = livestock,
        paddy = paddy,
        nonpaddy = nonpaddy,
        variables = vars,
    )
    return demand
end

"Initialize water demand for a water use sector"
function NonIrrigationDemand(nc, config, inds, dt; sector)
    demand_gross =
        ncread(
            nc,
            config,
            "vertical.demand.$(sector).demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.demand.$(sector).demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)
    n = length(inds)
    returnflow_f = return_flow_fraction.(demand_gross, demand_net)

    demand = PrescibedDemand{Float}(; demand_gross = demand_gross, demand_net = demand_net)
    vars = NonIrrigationDemandVariables{Float}(;
        returnflow_fraction = returnflow_f,
        returnflow = fill(Float(0), n),
    )
    non_irrigation_demand = NonIrrigationDemand{Float}(; demand = demand, variables = vars)

    return non_irrigation_demand
end

"Initialize paddy (rice) fields for water demand and irrigation computations"
function Paddy(nc, config, inds, dt)
    h_min = ncread(
        nc,
        config,
        "vertical.demand.paddy.parameters.h_min";
        sel = inds,
        defaults = 20.0,
        type = Float,
    )
    h_opt = ncread(
        nc,
        config,
        "vertical.demand.paddy.parameters.h_opt";
        sel = inds,
        defaults = 50.0,
        type = Float,
    )
    h_max = ncread(
        nc,
        config,
        "vertical.demand.paddy.parameters.h_max";
        sel = inds,
        defaults = 80.0,
        type = Float,
    )
    efficiency = ncread(
        nc,
        config,
        "vertical.demand.paddy.parameters.irrigation_efficiency";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    areas = ncread(
        nc,
        config,
        "vertical.demand.paddy.parameters.irrigation_areas";
        sel = inds,
        optional = false,
        type = Bool,
    )
    irrigation_trigger = ncread(
        nc,
        config,
        "vertical.demand.paddy.parameters.irrigation_trigger";
        sel = inds,
        optional = false,
        type = Bool,
    )
    max_irri_rate =
        ncread(
            nc,
            config,
            "vertical.demand.paddy.parameters.maximum_irrigation_rate";
            sel = inds,
            defaults = 25.0,
            type = Float,
        ) .* (dt / basetimestep)
    n = length(inds)
    params = PaddyParameters{Float}(;
        irrigation_efficiency = efficiency,
        maximum_irrigation_rate = max_irri_rate,
        irrigation_trigger = irrigation_trigger,
        h_min = h_min,
        h_max = h_max,
        h_opt = h_opt,
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

"Initialize crop (non paddy) fields for water demand and irrigation computations"
function NonPaddy(nc, config, inds, dt)
    efficiency = ncread(
        nc,
        config,
        "vertical.demand.nonpaddy.parameters.irrigation_efficiency";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    areas = ncread(
        nc,
        config,
        "vertical.demand.nonpaddy.parameters.irrigation_areas";
        sel = inds,
        defaults = 1,
        optional = false,
        type = Int,
    )
    irrigation_trigger = ncread(
        nc,
        config,
        "vertical.demand.nonpaddy.parameters.irrigation_trigger";
        sel = inds,
        defaults = 1,
        optional = false,
        type = Bool,
    )
    max_irri_rate =
        ncread(
            nc,
            config,
            "vertical.demand.nonpaddy.parameters.maximum_irrigation_rate";
            sel = inds,
            defaults = 25.0,
            type = Float,
        ) .* (dt / basetimestep)

    params = NonPaddyParameters{Float}(;
        maximum_irrigation_rate = max_irri_rate,
        irrigation_efficiency = efficiency,
        irrigation_areas = areas,
        irrigation_trigger = irrigation_trigger,
    )
    vars = NonPaddyVariables{Float}(; demand_gross = fill(mv, length(inds)))

    nonpaddy = NonPaddy{Float}(; variables = vars, parameters = params)

    return nonpaddy
end

"Initialize water allocation for the river domain"
function AllocationRiver(n)
    vars = AllocationRiverVariables{Float}(;
        act_surfacewater_abst = zeros(Float, n),
        act_surfacewater_abst_vol = zeros(Float, n),
        available_surfacewater = zeros(Float, n),
        nonirri_returnflow = zeros(Float, n),
    )
    allocation = AllocationRiver{Float}(; variables = vars)
    return allocation
end

"Initialize water allocation for the land domain (`vertical`)"
function AllocationLand(nc, config, inds)
    frac_sw_used = ncread(
        nc,
        config,
        "vertical.allocation.parameters.frac_sw_used";
        sel = inds,
        defaults = 1,
        type = Float,
    )
    areas = ncread(
        nc,
        config,
        "vertical.allocation.parameters.areas";
        sel = inds,
        defaults = 1,
        type = Int,
    )

    n = length(inds)

    params = AllocationLandParameters(; areas = areas, frac_sw_used = frac_sw_used)
    vars = AllocationLandVariables(;
        surfacewater_alloc = zeros(Float, n),
        act_groundwater_abst = zeros(Float, n),
        act_groundwater_abst_vol = zeros(Float, n),
        available_groundwater = zeros(Float, n),
        groundwater_alloc = zeros(Float, n),
        irri_alloc = zeros(Float, n),
        nonirri_alloc = zeros(Float, n),
        total_alloc = zeros(Float, n),
        nonirri_returnflow = zeros(Float, n),
    )
    allocation = AllocationLand(; parameters = params, variables = vars)
    return allocation
end

"Update returnflow fraction"
function return_flow_fraction!(model::NonIrrigationDemand)
    (; returnflow_fraction) = model.variables
    (; demand_gross, demand_net) = model.demand
    @. returnflow_fraction = return_flow_fraction(demand_gross, demand_net)
end

# return zero (gross demand) if non-irrigation sector is not defined
return_flow_fraction!(model::NoNonIrrigationDemand) = nothing

"Update water allocation for river and land domains based on local surface water (river) availability."
function surface_water_allocation_local(land_allocation, demand, river, network, dt)
    (; surfacewater_alloc) = land_allocation.variables
    (; surfacewater_demand) = demand.variables
    (; act_surfacewater_abst_vol, act_surfacewater_abst, available_surfacewater) =
        river.allocation.variables
    # maps from the land domain to the internal river domain (linear index), excluding water bodies
    index_river = network.land.index_river_wb
    for i in eachindex(surfacewater_demand)
        if index_river[i] > 0.0
            # the available volume is limited by a fixed scaling factor of 0.8 to prevent
            # rivers completely drying out. check for abstraction through inflow (external
            # negative inflow) and adjust available volume.
            if river.inflow[index_river[i]] < 0.0
                inflow = river.inflow[index_river[i]] * dt
                available_volume = max(river.volume[index_river[i]] * 0.80 + inflow, 0.0)
            else
                available_volume = river.volume[index_river[i]] * 0.80
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
end

"Update water allocation for river and land domains based on surface water (river) availability for allocation areas."
function surface_water_allocation_area(land_allocation, demand, river, network)
    inds_river = network.river.indices_allocation_areas
    inds_land = network.land.indices_allocation_areas
    res_index = network.river.reservoir_index
    lake_index = network.river.lake_index

    (; available_surfacewater, act_surfacewater_abst_vol, act_surfacewater_abst) =
        river.allocation.variables
    (; surfacewater_alloc) = land_allocation.variables
    (; surfacewater_demand) = demand.variables

    m = length(inds_river)
    # loop over allocation areas
    for i in 1:m
        # surface water demand (allocation area)
        sw_demand_vol = 0.0
        for j in inds_land[i]
            sw_demand_vol += surfacewater_demand[j] * 0.001 * network.land.area[j]
        end
        # surface water availability (allocation area)
        sw_available = 0.0
        for j in inds_river[i]
            if res_index[j] > 0
                # for reservoir locations use reservoir volume
                k = res_index[j]
                available_surfacewater[j] = river.reservoir.volume[k] * 0.98 # limit available reservoir volume
                sw_available += available_surfacewater[j]
            elseif lake_index[j] > 0
                # for lake locations use lake volume
                k = lake_index[j]
                available_surfacewater[j] = river.lake.storage[k] * 0.98 # limit available lake volume
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
                (act_surfacewater_abst_vol[j] / network.river.area[j]) * 1000.0
        end

        # water allocated to each land cell.
        for j in inds_land[i]
            surfacewater_alloc[j] += frac_allocate_sw * surfacewater_demand[j]
        end
    end
end

"Update water allocation for land domain based on local groundwater availability."
function groundwater_allocation_local(land_allocation, demand, groundwater_volume, network)
    (;
        surfacewater_alloc,
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = land_allocation.variables
    (; groundwater_demand, total_gross_demand) = demand.variables

    for i in eachindex(groundwater_demand)
        # groundwater demand based on allocation from surface water.
        groundwater_demand[i] = max(total_gross_demand[i] - surfacewater_alloc[i], 0.0)
        # land index excluding water bodies
        if network.index_wb[i]
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
end

"Update water allocation for land domain based on groundwater availability for allocation areas."
function groundwater_allocation_area(land_allocation, demand, network)
    inds_river = network.river.indices_allocation_areas
    inds_land = network.land.indices_allocation_areas
    (;
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = land_allocation.variables

    (; groundwater_demand) = demand.variables

    m = length(inds_river)
    # loop over allocation areas
    for i in 1:m
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
end

"Return and update non-irrigation sector (domestic, livestock, industry) return flow"
function return_flow(non_irri::NonIrrigationDemand, nonirri_demand_gross, nonirri_alloc)
    for i in eachindex(non_irri.variables.returnflow)
        frac = bounded_divide(non_irri.demand.demand_gross[i], nonirri_demand_gross[i])
        allocate = frac * nonirri_alloc[i]
        non_irri.variables.returnflow[i] =
            non_irri.variables.returnflow_fraction[i] * allocate
    end
    return non_irri.variables.returnflow
end

# return zero (return flow) if non-irrigation sector is not defined
return_flow(non_irri::NoNonIrrigationDemand, nonirri_demand_gross, nonirri_alloc) = 0.0

groundwater_volume(model::LateralSSF) = model.volume
groundwater_volume(model) = model.flow.aquifer.volume

"""
    update_water_allocation!(allocation, demand, lateral, network, dt)

Update water `allocation` (land domain) and water allocation of river domain (part of
`lateral`) based on `demand` for a single timestep. First, surface water abstraction is
computed to satisfy local water demand (non-irrigation and irrigation), and then updated
(including lakes and reservoirs) to satisfy the remaining water demand for allocation areas.
Then groundwater abstraction is computed to satisfy the remaining local water demand, and
then updated to satisfy the remaining water demand for allocation areas. Finally,
non-irrigation return flows are updated.
"""
function update_water_allocation!(land_allocation, demand::Demand, lateral, network, dt)
    river = lateral.river
    index_river = network.land.index_river_wb
    res_index_f = network.river.reservoir_index_f
    lake_index_f = network.river.lake_index_f
    (;
        groundwater_alloc,
        surfacewater_alloc,
        act_groundwater_abst,
        act_groundwater_abst_vol,
        total_alloc,
        irri_alloc,
        nonirri_alloc,
        nonirri_returnflow,
    ) = land_allocation.variables

    (; surfacewater_demand, nonirri_demand_gross, irri_demand_gross, total_gross_demand) =
        demand.variables

    (; frac_sw_used) = land_allocation.parameters
    (; act_surfacewater_abst, act_surfacewater_abst_vol) = river.allocation.variables

    surfacewater_alloc .= 0.0
    act_surfacewater_abst .= 0.0
    act_surfacewater_abst_vol .= 0.0
    # total surface water demand for each land cell
    @. surfacewater_demand =
        frac_sw_used * nonirri_demand_gross + frac_sw_used * irri_demand_gross

    # local surface water demand and allocation (river, excluding reservoirs and lakes)
    surface_water_allocation_local(land_allocation, demand, river, network, dt)
    # surface water demand and allocation for areas
    surface_water_allocation_area(land_allocation, demand, river, network)

    @. river.abstraction = act_surfacewater_abst_vol / dt

    # for reservoir and lake locations set river abstraction at zero and abstract volume
    # from reservoir and lake, including an update of lake waterlevel
    if !isnothing(river.reservoir)
        @. river.abstraction[res_index_f] = 0.0
        @. river.reservoir.volume -= act_surfacewater_abst_vol[res_index_f]
    elseif !isnothing(river.lake)
        @. river.abstraction[lake_index_f] = 0.0
        lakes = river.lake
        @. lakes.storage -= act_surfacewater_abst_vol[lake_index_f]
        @. lakes.waterlevel =
            waterlevel(lakes.storfunc, lakes.area, lakes.storage, lakes.sh)
    end

    groundwater_alloc .= 0.0
    act_groundwater_abst_vol .= 0.0
    act_groundwater_abst .= 0.0
    # local groundwater demand and allocation
    groundwater_allocation_local(
        land_allocation,
        demand,
        groundwater_volume(lateral.subsurface),
        network.land,
    )
    # groundwater demand and allocation for areas
    groundwater_allocation_area(land_allocation, demand, network)

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
update_water_allocation!(allocation, demand::NoDemand, lateral, network, dt) = nothing

function update_demand_gross!(demand::Demand)
    (; nonpaddy, paddy, domestic, industry, livestock) = demand
    (; irri_demand_gross, nonirri_demand_gross, total_gross_demand) = demand.variables
    # update gross water demands
    industry_dem = get_demand_gross(industry)
    domestic_dem = get_demand_gross(domestic)
    livestock_dem = get_demand_gross(livestock)
    nonpaddy_dem_gross = get_demand_gross(nonpaddy)
    paddy_dem_gross = get_demand_gross(paddy)
    @. irri_demand_gross = nonpaddy_dem_gross + paddy_dem_gross
    @. nonirri_demand_gross = industry_dem + domestic_dem + livestock_dem
    @. total_gross_demand =
        nonpaddy_dem_gross + paddy_dem_gross + industry_dem + domestic_dem + livestock_dem
end
update_demand_gross!(demand::NoDemand) = nothing

"""
    update_water_demand!(land_allocation, demand, soil::AbstractSoilModel)

Update water `demand` containing sectors `industry`, `domestic` and `livestock`, and `paddy`
rice fields and `nonpaddy` (other crop) fields.

Gross water demand for irrigation `irri_demand_gross` and non-irrigation
`nonirri_demand_gross`, and total gross water demand `total_gross_demand` are updated as
part of water allocation (`land_allocation`).
"""
function update_water_demand!(demand::Demand, soil)
    (; nonpaddy, paddy, domestic, industry, livestock) = demand

    return_flow_fraction!(industry)
    return_flow_fraction!(domestic)
    return_flow_fraction!(livestock)

    update_demand_gross!(nonpaddy, soil)
    update_demand_gross!(paddy)
    update_demand_gross!(demand)
end
update_water_demand!(demand::NoDemand, soil) = nothing