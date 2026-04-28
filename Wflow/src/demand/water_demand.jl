abstract type AbstractIrrigationModel end

struct NoIrrigationPaddyModel <: AbstractIrrigationModel
    n_land_cells::Int
end
struct NoIrrigationNonPaddyModel <: AbstractIrrigationModel
    n_land_cells::Int
end
struct NoNonIrrigationDemandModel <: AbstractDemandModel
    n_land_cells::Int
end
struct NoAllocationLandModel <: AbstractAllocationModel
    n_land_cells::Int
end
struct NoAllocationRiverModel <: AbstractAllocationModel
    n_land_cells::Int
end

"Struct to store non-irrigation water demand variables"
@with_kw struct NonIrrigationDemandVariables
    returnflow::Vector{Float64}               # return flow [mm Δt⁻¹]
    returnflow_fraction::Vector{Float64}      # return flow fraction [-]
end

"Struct to store prescribed water demand variables"
@with_kw struct PrescribedDemand
    demand_gross::Vector{Float64}     # gross water demand [mm Δt⁻¹]
    demand_net::Vector{Float64}       # net water demand [mm Δt⁻¹]
end

"Non-irrigation water demand model"
@with_kw struct NonIrrigationDemandModel <: AbstractDemandModel
    demand::PrescribedDemand
    variables::NonIrrigationDemandVariables
end

# wrapper methods
get_demand_gross(demand_model::NonIrrigationDemandModel) = demand_model.demand.demand_gross
get_demand_gross(demand_model::NoNonIrrigationDemandModel) =
    Zeros(demand_model.n_land_cells)

"Initialize non-irrigation water demand model for a water use `sector`"
function NonIrrigationDemandModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
    dt::Second,
    sector::AbstractString,
)
    demand_gross =
        ncread(
            dataset,
            config,
            "$(sector)__gross_water_demand_volume_flux",
            LandHydrologySBM;
            sel = land_indices_2d,
        ) .* (dt / BASETIMESTEP)
    demand_net =
        ncread(
            dataset,
            config,
            "$(sector)__net_water_demand_volume_flux",
            LandHydrologySBM;
            sel = land_indices_2d,
        ) .* (dt / BASETIMESTEP)
    n_land_cells = length(land_indices_2d)
    returnflow_f = return_flow_fraction.(demand_gross, demand_net)

    demand = PrescribedDemand(; demand_gross, demand_net)
    vars = NonIrrigationDemandVariables(;
        returnflow_fraction = returnflow_f,
        returnflow = zeros(n_land_cells),
    )
    non_irrigation_demand = NonIrrigationDemandModel(; demand, variables = vars)

    return non_irrigation_demand
end

"Struct to store non-paddy irrigation model variables"
@with_kw struct NonPaddyVariables
    n_land_cells::Int
    demand_gross::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)     # irrigation gross demand [mm Δt⁻¹]
end

"Struct to store non-paddy irrigation model parameters"
@with_kw struct NonPaddyParameters
    irrigation_efficiency::Vector{Float64}        # irrigation efficiency [-]
    maximum_irrigation_rate::Vector{Float64}      # maximum irrigation rate [mm Δt⁻¹]
    irrigation_areas::Vector{Bool}          # irrigation areas [-]
    irrigation_trigger::Vector{Bool}        # irrigation on or off [-]
end

"Non-paddy (other crops than flooded rice) irrigation model"
@with_kw struct NonPaddyModel <: AbstractIrrigationModel
    parameters::NonPaddyParameters
    variables::NonPaddyVariables
end

"Initialize non-paddy irrigation model"
function NonPaddyModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
    dt::Second,
)
    efficiency = ncread(
        dataset,
        config,
        "irrigated_non_paddy__irrigation_efficiency",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    areas = ncread(
        dataset,
        config,
        "irrigated_non_paddy_area__count",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    irrigation_trigger = ncread(
        dataset,
        config,
        "irrigated_non_paddy__irrigation_trigger_flag",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    max_irri_rate =
        ncread(
            dataset,
            config,
            "irrigated_non_paddy__max_irrigation_rate",
            LandHydrologySBM;
            sel = land_indices_2d,
        ) .* (dt / BASETIMESTEP)

    parameters = NonPaddyParameters(;
        maximum_irrigation_rate = max_irri_rate,
        irrigation_efficiency = efficiency,
        irrigation_areas = areas,
        irrigation_trigger,
    )
    n_land_cells = length(land_indices_2d)
    variables = NonPaddyVariables(; n_land_cells)

    nonpaddy_model = NonPaddyModel(; variables, parameters)

    return nonpaddy_model
end

# wrapper methods
get_demand_gross(nonpaddy_model::NonPaddyModel) = nonpaddy_model.variables.demand_gross
get_demand_gross(irrigation_model::NoIrrigationNonPaddyModel) =
    Zeros(irrigation_model.n_land_cells)

"""
    update_demand_gross!(nonpaddy_model::NonPaddyModel, soil_model::SbmSoilModel)

Update gross water demand `demand_gross` of the non-paddy irrigation model for a single
timestep.

The gross water demand is based on irrigation that is applied when the `irrigation_trigger`
is `true` (`on`) and when water depletion exceeds the readily available water in the root
zone of the SBM soil model. Irrigation brings the root zone back to field capacity, limited
by the infiltration capacity, taking into account limited irrigation efficiency and limited
by a maximum irrigation rate.
"""
function update_demand_gross!(nonpaddy_model::NonPaddyModel, soil_model::SbmSoilModel)
    (; parameters, variables) = nonpaddy_model
    (; irrigation_areas, irrigation_trigger, maximum_irrigation_rate) = parameters
    (; demand_gross, n_land_cells) = variables
    (; n_unsatlayers) = soil_model.variables

    for land_cell_idx in 1:n_land_cells
        if irrigation_areas[land_cell_idx] && irrigation_trigger[land_cell_idx]
            irri_dem_gross = 0.0
            for soil_layer_idx in 1:n_unsatlayers[land_cell_idx]
                depletion, readily_available_water =
                    water_demand_root_zone(soil_model, land_cell_idx, soil_layer_idx)

                # check if maximum irrigation rate has been applied at the previous time step.
                max_irri_rate_applied =
                    demand_gross[land_cell_idx] == maximum_irrigation_rate[land_cell_idx]
                if depletion >= readily_available_water # start irrigation
                    irri_dem_gross += depletion
                    # add depletion to irrigation gross demand when the maximum irrigation rate has been
                    # applied at the previous time step (to get volumetric water content at field capacity)
                elseif depletion > 0.0 && max_irri_rate_applied # continue irrigation
                    irri_dem_gross += depletion
                end
            end
            demand_gross[land_cell_idx] = compute_demand_gross(
                nonpaddy_model,
                soil_model,
                irri_dem_gross,
                land_cell_idx,
            )
        else
            demand_gross[land_cell_idx] = 0.0
        end
    end
    return nothing
end

update_demand_gross!(nonpaddy_model::NoIrrigationNonPaddyModel, soil_model::SbmSoilModel) =
    nothing

"Compute water demand only for root zone through root fraction per layer"
function water_demand_root_zone(
    soil_model::SbmSoilModel,
    land_cell_idx::Int,
    soil_layer_idx::Int,
)
    (; sumlayers, hb, theta_s, theta_r, theta_fc, c) = soil_model.parameters
    (; ustorelayerthickness, ustorelayerdepth, h3) = soil_model.variables

    rootingdepth = get_rootingdepth(soil_model)

    rootfrac = min(
        1.0,
        (
            max(
                0.0,
                rootingdepth[land_cell_idx] - sumlayers[land_cell_idx][soil_layer_idx],
            ) / ustorelayerthickness[land_cell_idx][soil_layer_idx]
        ),
    )
    # vwc_h3 can be precalculated.
    vwc_h3 = vwc_brooks_corey(
        h3[land_cell_idx],
        hb[land_cell_idx],
        theta_s[land_cell_idx],
        theta_r[land_cell_idx],
        c[land_cell_idx][soil_layer_idx],
    )
    depletion =
        (theta_fc[land_cell_idx] * ustorelayerthickness[land_cell_idx][soil_layer_idx]) - (
            ustorelayerdepth[land_cell_idx][soil_layer_idx] +
            theta_r[land_cell_idx] * ustorelayerthickness[land_cell_idx][soil_layer_idx]
        )
    depletion *= rootfrac
    readily_available_water =
        (theta_fc[land_cell_idx] - vwc_h3) *
        ustorelayerthickness[land_cell_idx][soil_layer_idx]
    readily_available_water *= rootfrac

    return depletion, readily_available_water
end

function compute_demand_gross(
    nonpaddy_model::NonPaddyModel,
    soil_model::SbmSoilModel,
    irri_dem_gross::Float64,
    land_cell_idx::Int,
)
    (; pathfrac, infiltcapsoil) = soil_model.parameters
    (; f_infiltration_reduction) = soil_model.variables
    (; irrigation_efficiency, maximum_irrigation_rate) = nonpaddy_model.parameters

    infiltration_capacity =
        f_infiltration_reduction[land_cell_idx] *
        (1.0 - pathfrac[land_cell_idx]) *
        infiltcapsoil[land_cell_idx]
    irri_dem_gross = min(irri_dem_gross, infiltration_capacity)
    irri_dem_gross /= irrigation_efficiency[land_cell_idx]
    # limit irrigation demand to the maximum irrigation rate
    irri_dem_gross = min(irri_dem_gross, maximum_irrigation_rate[land_cell_idx])

    return irri_dem_gross
end

"Struct to store paddy irrigation model variables"
@with_kw struct PaddyVariables
    n_land_cells::Int
    demand_gross::Vector{Float64} = fill(MISSING_VALUE, n_land_cells) # irrigation gross demand [mm Δt⁻¹]
    h::Vector{Float64} = zeros(n_land_cells)                          # actual water depth in rice field [mm]
    evaporation::Vector{Float64} = zeros(n_land_cells)                # evaporation rate [mm Δt⁻¹]
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

"PaddyModel (flooded rice) irrigation model"
@with_kw struct PaddyModel <: AbstractIrrigationModel
    parameters::PaddyParameters
    variables::PaddyVariables
end

"Initialize paddy irrigation model"
function PaddyModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
    dt::Second,
)
    h_min = ncread(
        dataset,
        config,
        "irrigated_paddy__min_depth",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    h_opt = ncread(
        dataset,
        config,
        "irrigated_paddy__optimal_depth",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    h_max = ncread(
        dataset,
        config,
        "irrigated_paddy__max_depth",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    efficiency = ncread(
        dataset,
        config,
        "irrigated_paddy__irrigation_efficiency",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    areas = ncread(
        dataset,
        config,
        "irrigated_paddy_area__count",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    irrigation_trigger = ncread(
        dataset,
        config,
        "irrigated_paddy__irrigation_trigger_flag",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    max_irri_rate =
        ncread(
            dataset,
            config,
            "irrigated_paddy__max_irrigation_rate",
            LandHydrologySBM;
            sel = land_indices_2d,
        ) .* (dt / BASETIMESTEP)
    n_land_cells = length(land_indices_2d)
    parameters = PaddyParameters(;
        irrigation_efficiency = efficiency,
        maximum_irrigation_rate = max_irri_rate,
        irrigation_trigger,
        h_min,
        h_max,
        h_opt,
        irrigation_areas = areas,
    )
    variables = PaddyVariables(; n_land_cells)
    paddy = PaddyModel(; parameters, variables)
    return paddy
end

# wrapper methods
get_water_depth(paddy_model::PaddyModel) = paddy_model.variables.h
get_water_depth(paddy_model::NoIrrigationPaddyModel) = Zeros(paddy_model.n_land_cells)
get_demand_gross(paddy_model::PaddyModel) = paddy_model.variables.demand_gross
get_demand_gross(paddy_model::NoIrrigationPaddyModel) = Zeros(paddy_model.n_land_cells)

"""
    evaporation!(paddy_model::PaddyModel, potential_evaporation)

Update `evaporation` and the water depth `h` of the paddy irrigation model for a single
timestep.
"""
function evaporation!(paddy_model::PaddyModel, potential_evaporation)
    (; n_land_cells) = paddy_model.variables

    for land_cell_idx in 1:n_land_cells
        if paddy_model.parameters.irrigation_areas[land_cell_idx]
            evaporation = min(
                paddy_model.variables.h[land_cell_idx],
                potential_evaporation[land_cell_idx],
            )
            paddy_model.variables.h[land_cell_idx] -= evaporation
            paddy_model.variables.evaporation[land_cell_idx] = evaporation
        end
    end
    return nothing
end
evaporation!(::NoIrrigationPaddyModel, ::Any) = nothing

# wrapper methods
get_evaporation(paddy_model::NoIrrigationPaddyModel) = Zeros(paddy_model.n_land_cells)
get_evaporation(paddy_model::PaddyModel) = paddy_model.variables.evaporation

"""
    update_runoff!(paddy_model::PaddyModel, runoff)

Update `runoff` based on the water depth `h_max` (paddy field starts spilling), and update
the water depth `h` of the paddy irrigation model for a single timestep.
"""
function update_runoff!(paddy_model::PaddyModel, runoff)
    (; n_land_cells) = paddy_model.variables
    for land_cell_idx in 1:n_land_cells
        if paddy_model.parameters.irrigation_areas[land_cell_idx]
            paddy_runoff = max(
                runoff[land_cell_idx] - paddy_model.parameters.h_max[land_cell_idx],
                0.0,
            )
            paddy_model.variables.h[land_cell_idx] = runoff[land_cell_idx] - paddy_runoff
            runoff[land_cell_idx] = paddy_runoff
        end
    end
    return nothing
end
update_runoff!(::NoIrrigationPaddyModel, ::Any) = nothing

"""
    update_demand_gross!(paddy_model::PaddyModel)

Update gross water demand `demand_gross` of the paddy irrigation model for a single
timestep.

The gross water demand is based on irrigation that is applied when the `irrigation_trigger`
is `true` (`on`) and when the paddy water depth `h` reaches below the minimum water depth
`h_min`. Irrigation is the amount required to reach the optimal paddy water depth `h_opt`,
taking into account limited irrigation efficiency and limited by a maximum irrigation rate.
"""
function update_demand_gross!(paddy_model::PaddyModel)
    (; demand_gross, n_land_cells) = paddy_model.variables
    (;
        irrigation_areas,
        irrigation_trigger,
        irrigation_efficiency,
        maximum_irrigation_rate,
    ) = paddy_model.parameters

    for land_cell_idx in 1:n_land_cells
        if irrigation_areas[land_cell_idx] && irrigation_trigger[land_cell_idx]
            irr_depth_paddy = compute_irrigation_depth(paddy_model, land_cell_idx)

            irri_dem_gross = irr_depth_paddy / irrigation_efficiency[land_cell_idx]
            # limit irrigation demand to the maximum irrigation rate
            irri_dem_gross = min(irri_dem_gross, maximum_irrigation_rate[land_cell_idx])
            demand_gross[land_cell_idx] = irri_dem_gross
        else
            demand_gross[land_cell_idx] = 0.0
        end
    end
end

function compute_irrigation_depth(paddy_model::PaddyModel, land_cell_idx::Int)
    (; maximum_irrigation_rate, h_min, h_opt) = paddy_model.parameters
    (; demand_gross, h) = paddy_model.variables

    # check if maximum irrigation rate has been applied at the previous time step.
    max_irri_rate_applied =
        demand_gross[land_cell_idx] == maximum_irrigation_rate[land_cell_idx]
    # start irrigation
    irr_depth_paddy = if h[land_cell_idx] < h_min[land_cell_idx]
        h_opt[land_cell_idx] - h[land_cell_idx]
    elseif h[land_cell_idx] < h_opt[land_cell_idx] && max_irri_rate_applied # continue irrigation
        h_opt[land_cell_idx] - h[land_cell_idx]
    else
        0.0
    end

    return irr_depth_paddy
end

update_demand_gross!(::NoIrrigationPaddyModel) = nothing

"Struct to store water demand model variables"
@with_kw struct DemandVariables
    n_land_cells::Int
    irri_demand_gross::Vector{Float64} = zeros(n_land_cells)        # irrigation gross demand [mm Δt⁻¹]
    nonirri_demand_gross::Vector{Float64} = zeros(n_land_cells)     # non-irrigation gross demand [mm Δt⁻¹]
    total_gross_demand::Vector{Float64} = zeros(n_land_cells)       # total gross demand [mm Δt⁻¹]
    surfacewater_demand::Vector{Float64} = zeros(n_land_cells)      # demand from surface water [mm Δt⁻¹]
    groundwater_demand::Vector{Float64} = zeros(n_land_cells)       # demand from groundwater [mm Δt⁻¹]
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
    n_land_cells::Int
    domestic::NoNonIrrigationDemandModel = NoNonIrrigationDemandModel(n_land_cells)
    industry::NoNonIrrigationDemandModel = NoNonIrrigationDemandModel(n_land_cells)
    livestock::NoNonIrrigationDemandModel = NoNonIrrigationDemandModel(n_land_cells)
    paddy::NoIrrigationPaddyModel = NoIrrigationPaddyModel(n_land_cells)
    nonpaddy::NoIrrigationNonPaddyModel = NoIrrigationNonPaddyModel(n_land_cells)
end

"Initialize water demand model"
function DemandModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
    dt::Second,
)
    n_land_cells = length(land_indices_2d)
    demand(
        name;
        constr = NonIrrigationDemandModel,
        constr_triv = NoNonIrrigationDemandModel,
    ) =
        if getfield(config.model.water_demand, Symbol("$(name)__flag"))::Bool
            if constr == NonIrrigationDemandModel
                constr(dataset, config, land_indices_2d, dt, name)
            else
                constr(dataset, config, land_indices_2d, dt)
            end
        else
            constr_triv(n_land_cells)
        end

    domestic = demand("domestic")
    industry = demand("industry")
    livestock = demand("livestock")
    paddy = demand("paddy"; constr = PaddyModel, constr_triv = NoIrrigationPaddyModel)
    nonpaddy =
        demand("nonpaddy"; constr = NonPaddyModel, constr_triv = NoIrrigationNonPaddyModel)

    variables = DemandVariables(; n_land_cells)
    demand_model = DemandModel(; domestic, industry, livestock, paddy, nonpaddy, variables)
    return demand_model
end

"Struct to store river allocation model variables"
@with_kw struct AllocationRiverVariables
    n_river_cells::Int
    act_surfacewater_abst::Vector{Float64} = zeros(n_river_cells)        # actual surface water abstraction [mm Δt⁻¹]
    act_surfacewater_abst_vol::Vector{Float64} = zeros(n_river_cells)    # actual surface water abstraction [m³ Δt⁻¹]
    available_surfacewater::Vector{Float64} = zeros(n_river_cells)       # available surface water [m³]
    nonirri_returnflow::Vector{Float64} = zeros(n_river_cells)           # return flow from non irrigation [mm Δt⁻¹]
end

"River allocation model"
@with_kw struct AllocationRiverModel <: AbstractAllocationModel
    n_river_cells::Int
    variables::AllocationRiverVariables = AllocationRiverVariables(; n_river_cells)
end

get_nonirrigation_returnflow(allocation_model::AllocationRiverModel) =
    allocation_model.variables.nonirri_returnflow
get_nonirrigation_returnflow(allocation_model::NoAllocationRiverModel) =
    Zeros(allocation_model.n_land_cells)

"Struct to store land allocation allocation model parameters"
@with_kw struct AllocationLandParameters
    frac_sw_used::Vector{Float64}     # fraction surface water used [-]
    areas::Vector{Int}          # allocation areas [-]
end

"Struct to store land allocation model variables"
@with_kw struct AllocationLandVariables
    n_land_cells::Int
    surfacewater_alloc::Vector{Float64} = zeros(n_land_cells)           # allocation from surface water [mm Δt⁻¹]
    act_groundwater_abst::Vector{Float64} = zeros(n_land_cells)         # actual groundwater abstraction [mm Δt⁻¹]
    act_groundwater_abst_vol::Vector{Float64} = zeros(n_land_cells)     # actual groundwater abstraction [m³ Δt⁻¹]
    available_groundwater::Vector{Float64} = zeros(n_land_cells)        # available groundwater [m³]
    groundwater_alloc::Vector{Float64} = zeros(n_land_cells)            # allocation from groundwater [mm Δt⁻¹]
    irri_alloc::Vector{Float64} = zeros(n_land_cells)                   # allocated water for irrigation [mm Δt⁻¹]
    nonirri_alloc::Vector{Float64} = zeros(n_land_cells)                # allocated water for non-irrigation [mm Δt⁻¹]
    total_alloc::Vector{Float64} = zeros(n_land_cells)                  # total allocated water [mm Δt⁻¹]
    nonirri_returnflow::Vector{Float64} = zeros(n_land_cells)           # return flow from non irrigation [mm Δt⁻¹]
end

"Land allocation model"
@with_kw struct AllocationLandModel <: AbstractAllocationModel
    n_land_cells::Int
    parameters::AllocationLandParameters
    variables::AllocationLandVariables = AllocationLandVariables(; n_land_cells)
end

"Initialize water allocation for the land domain"
function AllocationLandModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
)
    frac_sw_used = ncread(
        dataset,
        config,
        "land_surface_water__withdrawal_fraction",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    areas = ncread(
        dataset,
        config,
        "land_water_allocation_area__count",
        LandHydrologySBM;
        sel = land_indices_2d,
    )

    n_land_cells = length(land_indices_2d)

    parameters = AllocationLandParameters(; areas, frac_sw_used)
    allocation = AllocationLandModel(; n_land_cells, parameters)
    return allocation
end

# wrapper methods
get_irrigation_allocated(allocation_model::AllocationLandModel) =
    allocation_model.variables.irri_alloc
get_irrigation_allocated(allocation_model::NoAllocationLandModel) =
    Zeros(allocation_model.n_land_cells)
get_nonirrigation_returnflow(allocation_model::AllocationLandModel) =
    allocation_model.variables.nonirri_returnflow
get_nonirrigation_returnflow(allocation_model::NoAllocationLandModel) =
    Zeros(allocation_model.n_land_cells)
get_groundwater_abstraction_flux(allocation_model::AllocationLandModel) =
    allocation_model.variables.act_groundwater_abst
get_groundwater_abstraction_flux(allocation_model::NoAllocationLandModel) =
    Zeros(allocation_model.n_land_cells)

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
    (; surfacewater_alloc) = allocation_model.variables
    (; surfacewater_demand) = demand_variables
    (; act_surfacewater_abst_vol, act_surfacewater_abst, available_surfacewater) =
        river_flow_model.allocation.variables
    (; external_inflow) = river_flow_model.boundary_conditions
    (; storage) = river_flow_model.variables
    (; area) = domain.parameters
    land_cell_river_indices = domain.network.land_cell_river_indices_excl_reservoir

    # maps from the land domain to the internal river domain (linear index), excluding reservoirs
    for (land_cell_idx, river_cell_idx) in enumerate(land_cell_river_indices)
        if river_cell_idx > 0
            # the available volume is limited by a fixed scaling factor of 0.8 to prevent
            # rivers completely drying out. check for abstraction through negative external
            # inflow first and adjust available volume.
            available_volume = storage[river_cell_idx] * 0.80
            if external_inflow[river_cell_idx] < 0.0
                max_river_abstraction =
                    min(-external_inflow[river_cell_idx] * dt, available_volume)
                available_volume = max(available_volume - max_river_abstraction, 0.0)
            end
            # satisfy surface water demand with available local river volume
            surfacewater_demand_vol =
                surfacewater_demand[land_cell_idx] * 0.001 * area[land_cell_idx]
            abstraction_vol = min(surfacewater_demand_vol, available_volume)
            act_surfacewater_abst_vol[river_cell_idx] = abstraction_vol
            # remaining available surface water and demand
            available_surfacewater[river_cell_idx] =
                max(available_volume - abstraction_vol, 0.0)
            abstraction = (abstraction_vol / area[land_cell_idx]) * 1000.0
            surfacewater_demand[land_cell_idx] =
                max(surfacewater_demand[land_cell_idx] - abstraction, 0.0)
            # update actual abstraction from river and surface water allocation (land cell)
            act_surfacewater_abst[river_cell_idx] = abstraction
            surfacewater_alloc[land_cell_idx] = abstraction
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

    (; available_surfacewater, act_surfacewater_abst_vol, act_surfacewater_abst) =
        river_flow_model.allocation.variables
    (; surfacewater_alloc) = allocation_model.variables
    (; surfacewater_demand) = demand_variables
    (; reservoir) = river_flow_model.boundary_conditions

    for land_cell_idx in eachindex(inds_river)
        # surface water_demand (allocation area)
        sw_demand_vol = mapreduce(
            land_cell_idx_other ->
                surfacewater_demand[land_cell_idx_other] * 1e-3 * area[land_cell_idx_other],
            +,
            inds_land[land_cell_idx],
        )

        sw_available = available_surface_water!(
            available_surfacewater,
            reservoir,
            inds_river[land_cell_idx],
            inds_reservoir,
            dt,
        )

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
        for land_cell_idx_other in inds_river[land_cell_idx]
            act_surfacewater_abst_vol[land_cell_idx_other] +=
                frac_abstract_sw * available_surfacewater[land_cell_idx_other]
            act_surfacewater_abst[land_cell_idx_other] =
                (
                    act_surfacewater_abst_vol[land_cell_idx_other] /
                    domain.river.parameters.cell_area[land_cell_idx_other]
                ) * 1000.0
        end

        # water allocated to each land cell.
        for land_cell_idx_other in inds_land[land_cell_idx]
            surfacewater_alloc[land_cell_idx_other] +=
                frac_allocate_sw * surfacewater_demand[land_cell_idx_other]
        end
    end
end

function available_surface_water!(
    available_surfacewater::Vector{Float64},
    reservoir_model,
    indices_river::Vector{Int},
    indices_reservoir::Vector{Int},
    dt::Float64,
)
    sw_available = 0.0
    for land_cell_idx in indices_river
        reservoir_idx = indices_reservoir[land_cell_idx]
        if reservoir_idx > 0
            # for reservoir locations use reservoir storage, check for abstraction
            # through external negative inflow first and adjust available volume.
            external_inflow =
                reservoir_model.boundary_conditions.external_inflow[reservoir_idx]
            available_volume = reservoir_model.variables.storage[reservoir_idx] * 0.98
            if external_inflow < 0.0
                if available_volume > -external_inflow * dt
                    available_volume += external_inflow * dt
                else
                    available_volume = 0.0
                end
            end
            available_surfacewater[land_cell_idx] = available_volume
            sw_available += available_volume
        else
            # river volume
            sw_available += available_surfacewater[land_cell_idx]
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
)
    (;
        surfacewater_alloc,
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = allocation_model.variables
    (; groundwater_demand, total_gross_demand, n_land_cells) = demand_variables
    (; area, reservoir_coverage) = parameters

    for land_cell_idx in 1:n_land_cells
        # groundwater demand based on allocation from surface water.
        groundwater_demand[land_cell_idx] =
            max(total_gross_demand[land_cell_idx] - surfacewater_alloc[land_cell_idx], 0.0)
        # excluding reservoirs
        if !reservoir_coverage[land_cell_idx]
            # satisfy groundwater demand with available local groundwater volume
            groundwater_demand_vol =
                groundwater_demand[land_cell_idx] * 0.001 * area[land_cell_idx]
            available_volume = groundwater_storage[land_cell_idx] * 0.75 # limit available groundwater volume
            abstraction_vol = min(groundwater_demand_vol, available_volume)
            act_groundwater_abst_vol[land_cell_idx] = abstraction_vol
            # remaining available groundwater and demand
            available_groundwater[land_cell_idx] =
                max(available_volume - abstraction_vol, 0.0)
            abstraction = (abstraction_vol / area[land_cell_idx]) * 1000.0
            groundwater_demand[land_cell_idx] =
                max(groundwater_demand[land_cell_idx] - abstraction, 0.0)
            # update actual abstraction from groundwater and groundwater allocation (land cell)
            act_groundwater_abst[land_cell_idx] = abstraction
            groundwater_alloc[land_cell_idx] = abstraction
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
)
    inds_river = domain.river.network.allocation_area_indices
    inds_land = domain.land.network.allocation_area_indices
    (;
        act_groundwater_abst_vol,
        available_groundwater,
        act_groundwater_abst,
        groundwater_alloc,
    ) = allocation_model.variables

    (; groundwater_demand) = demand_variables
    (; area) = domain.land.parameters

    # loop over allocation areas
    for alloc_area_idx in eachindex(inds_river)
        # groundwater demand and availability (allocation area)
        gw_demand_vol = 0.0
        gw_available = 0.0
        for land_cell_idx in inds_land[alloc_area_idx]
            gw_demand_vol += groundwater_demand[land_cell_idx] * 0.001 * area[land_cell_idx]
            gw_available += available_groundwater[land_cell_idx]
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
        for land_cell_idx in inds_land[alloc_area_idx]
            act_groundwater_abst_vol[land_cell_idx] +=
                frac_abstract_gw * available_groundwater[land_cell_idx]
            act_groundwater_abst[land_cell_idx] =
                1000.0 * (act_groundwater_abst_vol[land_cell_idx] / area[land_cell_idx])
            groundwater_alloc[land_cell_idx] +=
                frac_allocate_gw * groundwater_demand[land_cell_idx]
        end
    end
    return nothing
end

"Return and update non-irrigation sector (domestic, livestock, industry) return flow"
function return_flow(
    demand_model::NonIrrigationDemandModel,
    nonirri_demand_gross::Vector{Float64},
    nonirri_alloc::Vector{Float64},
)
    for land_cell_idx in eachindex(demand_model.variables.returnflow)
        frac = bounded_divide(
            demand_model.demand.demand_gross[land_cell_idx],
            nonirri_demand_gross[land_cell_idx],
        )
        allocate = frac * nonirri_alloc[land_cell_idx]
        demand_model.variables.returnflow[land_cell_idx] =
            demand_model.variables.returnflow_fraction[land_cell_idx] * allocate
    end
    return demand_model.variables.returnflow
end

# return zero (return flow) if non-irrigation sector is not defined
return_flow(
    demand_model::NoNonIrrigationDemandModel,
    nonirri_demand_gross::Vector{Float64},
    nonirri_alloc::Vector{Float64},
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
    inds_river = domain.land.network.land_cell_river_indices_excl_reservoir
    inds_reservoir = domain.reservoir.network.river_cell_indices_containing_reservoir
    (;
        groundwater_alloc,
        surfacewater_alloc,
        act_groundwater_abst,
        act_groundwater_abst_vol,
        total_alloc,
        irri_alloc,
        nonirri_alloc,
        nonirri_returnflow,
        n_land_cells,
    ) = allocation_model.variables

    (; surfacewater_demand, nonirri_demand_gross, irri_demand_gross, total_gross_demand) =
        demand_model.variables

    (; frac_sw_used) = allocation_model.parameters
    (; act_surfacewater_abst, act_surfacewater_abst_vol) = river.allocation.variables
    (; abstraction, reservoir) = river.boundary_conditions

    surfacewater_alloc .= 0.0
    act_surfacewater_abst .= 0.0
    act_surfacewater_abst_vol .= 0.0
    # total surface water demand for each land cell
    @. surfacewater_demand =
        frac_sw_used * nonirri_demand_gross + frac_sw_used * irri_demand_gross

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
        allocation_model,
        demand_model.variables,
        routing.subsurface_flow.variables.storage,
        domain.land.parameters,
    )
    # groundwater demand and allocation for areas
    groundwater_allocation_area!(allocation_model, demand_model.variables, domain)

    # irrigation allocation
    for land_cell_idx in 1:n_land_cells
        total_alloc[land_cell_idx] =
            groundwater_alloc[land_cell_idx] + surfacewater_alloc[land_cell_idx]
        frac_irri = bounded_divide(
            irri_demand_gross[land_cell_idx],
            total_gross_demand[land_cell_idx],
        )
        irri_alloc[land_cell_idx] = frac_irri * total_alloc[land_cell_idx]
        nonirri_alloc[land_cell_idx] =
            total_alloc[land_cell_idx] - irri_alloc[land_cell_idx]
    end

    # non-irrigation return flows
    returnflow_livestock =
        return_flow(demand_model.livestock, nonirri_demand_gross, nonirri_alloc)
    returnflow_domestic =
        return_flow(demand_model.domestic, nonirri_demand_gross, nonirri_alloc)
    returnflow_industry =
        return_flow(demand_model.industry, nonirri_demand_gross, nonirri_alloc)

    @. nonirri_returnflow = returnflow_livestock + returnflow_domestic + returnflow_industry

    for land_cell_idx in eachindex(nonirri_returnflow)
        if inds_river[land_cell_idx] > 0
            river_cell_idx = inds_river[land_cell_idx]
            river.allocation.variables.nonirri_returnflow[river_cell_idx] =
                nonirri_returnflow[land_cell_idx]
            nonirri_returnflow[land_cell_idx] = 0.0
        else
            nonirri_returnflow[land_cell_idx] = nonirri_returnflow[land_cell_idx]
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
    update_demand_gross!(demand_model::DemandModel)

Update total irrigation gross water demand `irri_demand_gross`, total non-irrigation gross
water demand `nonirri_demand_gross` and total gross water demand `total_gross_demand`.
"""
function update_demand_gross!(demand_model::DemandModel)
    (; nonpaddy, paddy, domestic, industry, livestock) = demand_model
    (; irri_demand_gross, nonirri_demand_gross, total_gross_demand) = demand_model.variables
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

update_demand_gross!(demand_model::NoDemandModel) = nothing

"""
    update_water_demand_model!(demand_model::DemandModel, soil_model::SbmSoilModel)

Update the return flow fraction `returnflow_fraction` of `industry`, `domestic` and
`livestock`, gross water demand `demand_gross` of `paddy` and `nonpaddy` models, and the
total gross water demand, total irrigation gross water demand and total non-irrigation gross
water demand as part of the water `demand` model.
"""
function update_water_demand_model!(demand_model::DemandModel, soil_model::SbmSoilModel)
    (; nonpaddy, paddy, domestic, industry, livestock) = demand_model

    return_flow_fraction!(industry)
    return_flow_fraction!(domestic)
    return_flow_fraction!(livestock)

    update_demand_gross!(nonpaddy, soil_model)
    update_demand_gross!(paddy)
    update_demand_gross!(demand_model)

    return nothing
end
update_water_demand_model!(demand_model::NoDemandModel, soil_model::SbmSoilModel) = nothing
