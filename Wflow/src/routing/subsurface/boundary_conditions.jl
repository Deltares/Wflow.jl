function check_flux(flux::Float64, subsurface_flow_model::GroundwaterFlowModel, index::Int)
    # Check if cell is dry
    if subsurface_flow_model.variables.head[index] <=
       subsurface_flow_model.parameters.bottom[index]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end

function check_flux(flux::Float64, subsurface_flow_model::LateralSSFModel, index::Int)
    # Check if cell is dry
    if subsurface_flow_model.variables.zi[index] >=
       subsurface_flow_model.parameters.soilthickness[index]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end

@with_data_lookup struct GwfRiverParameters
    "river_water__infiltration_conductance"
    infiltration_conductance::Vector{Float64} # [m² d⁻¹]
    "river_water__exfiltration_conductance"
    exfiltration_conductance::Vector{Float64} # [m² d⁻¹]
    "river_bottom__elevation"
    bottom::Vector{Float64} # [m]
end

@with_data_lookup struct GwfRiverVariables
    n::Int
    stage::Vector{Float64} = fill(MISSING_VALUE, n) # [m]
    storage::Vector{Float64} = fill(MISSING_VALUE, n) # [m³]
    flux::Vector{Float64} = fill(MISSING_VALUE, n)  # [m³ d⁻¹]
    "river_water__to_subsurface_volume_flow_rate"
    flux_av::Vector{Float64} = fill(MISSING_VALUE, n)  # [m³ d⁻¹]
end

@kwdef struct GwfRiverModel <: AbstractSubsurfaceFlowBC
    parameters::GwfRiverParameters
    variables::GwfRiverVariables
end

function GwfRiverModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup=DataLookup(),
)
    infiltration_conductance = ncread(
        dataset,
        config,
        "river_water__infiltration_conductance",
        Routing;
        sel=indices,
    )
    exfiltration_conductance = ncread(
        dataset,
        config,
        "river_water__exfiltration_conductance",
        Routing;
        sel=indices,
    )
    bottom = ncread(dataset, config, "river_bottom__elevation", Routing; sel=indices)

    parameters = GwfRiverParameters(
        data_lookup;
        infiltration_conductance,
        exfiltration_conductance,
        bottom,
    )
    n = length(indices)
    variables = GwfRiverVariables(data_lookup; n)
    river_model = GwfRiverModel(parameters, variables)
    return river_model
end

function flux!(
    gwf_river_model::GwfRiverModel,
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        head = subsurface_flow_model.variables.head[index]
        stage = gwf_river_model.variables.stage[i]
        if stage > head
            max_infiltration_flux = gwf_river_model.variables.storage[i] / dt
            cond = gwf_river_model.parameters.infiltration_conductance[i]
            delta_head = min(stage - gwf_river_model.parameters.bottom[i], stage - head)
            flux = min(cond * delta_head, max_infiltration_flux)
        else
            cond = gwf_river_model.parameters.exfiltration_conductance[i]
            delta_head = stage - head
            flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        end
        gwf_river_model.variables.flux[i] = flux
        subsurface_flow_model.variables.q_net_bnds[index] += flux
        gwf_river_model.variables.storage[i] -= dt * flux
        gwf_river_model.variables.flux_av[i] += dt * flux
    end
    return nothing
end

@with_data_lookup struct DrainageParameters
    "land_drain__elevation"
    elevation::Vector{Float64} # [m]
    "land_drain__conductance"
    conductance::Vector{Float64} # [m² d⁻¹]
end

@with_data_lookup struct DrainageVariables
    n::Int
    flux::Vector{Float64} = fill(MISSING_VALUE, n) # [m³ d⁻¹]
    "land_drain_water__to_subsurface_volume_flow_rate"
    flux_av::Vector{Float64} = fill(MISSING_VALUE, n) # [m³ d⁻¹]
end

@kwdef struct DrainageModel <: AbstractSubsurfaceFlowBC
    parameters::DrainageParameters
    variables::DrainageVariables
end

function DrainageModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup=DataLookup(),
)
    elevation = ncread(dataset, config, "land_drain__elevation", Routing; sel=indices)
    conductance = ncread(dataset, config, "land_drain__conductance", Routing; sel=indices)
    parameters = DrainageParameters(data_lookup; elevation, conductance)
    n = length(indices)
    variables = DrainageVariables(data_lookup; n)

    drainage_model = DrainageModel(parameters, variables)
    return drainage_model
end

function flux!(
    drainage_model::DrainageModel,
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        cond = drainage_model.parameters.conductance[i]
        delta_head = min(
            0,
            drainage_model.parameters.elevation[i] -
            subsurface_flow_model.variables.head[index],
        )
        flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        drainage_model.variables.flux[i] = flux
        drainage_model.variables.flux_av[i] += dt * flux
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@kwdef struct HeadBoundaryParameters
    conductance::Vector{Float64} # [m² d⁻¹]
end

@kwdef struct HeadBoundaryVariables
    head::Vector{Float64} # [m]
    flux::Vector{Float64} # [m³ d⁻¹]
    flux_av::Vector{Float64} # [m³ d⁻¹]
end

@kwdef struct HeadBoundary <: AbstractSubsurfaceFlowBC
    parameters::HeadBoundaryParameters
    variables::HeadBoundaryVariables
end

function flux!(
    headboundary::HeadBoundary,
    subsurface_flow_model::GroundwaterFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        cond = headboundary.parameters.conductance[i]
        delta_head =
            headboundary.variables.head[i] - subsurface_flow_model.variables.head[index]
        flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        headboundary.variables.flux[i] = flux
        headboundary.variables.flux_av[i] += dt * flux
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@with_data_lookup struct RechargeVariables
    n::Int
    rate::Vector{Float64} = fill(MISSING_VALUE, n) # [m d⁻¹]
    flux::Vector{Float64} = zeros(n) # [m³ d⁻¹]
    "subsurface_water_saturated_zone_top__net_recharge_volume_flow_rate"
    flux_av::Vector{Float64} = zeros(n) # [m³ d⁻¹]
end

@kwdef struct RechargeModel <: AbstractSubsurfaceFlowBC
    n::Int
    variables::RechargeVariables = RechargeVariables(; n)
end

function flux!(
    recharge_model::RechargeModel,
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        flux = check_flux(
            recharge_model.variables.rate[i] * subsurface_flow_model.parameters.area[index],
            subsurface_flow_model,
            index,
        )
        recharge_model.variables.flux[i] = flux
        recharge_model.variables.flux_av[i] += dt * flux
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@kwdef struct WellVariables
    volumetric_rate::Vector{Float64} # [m³ d⁻¹]
    flux::Vector{Float64} # [m³ d⁻¹]
    flux_av::Vector{Float64} # [m³ d⁻¹]
end

@kwdef struct WellModel <: AbstractSubsurfaceFlowBC
    variables::WellVariables
end

function flux!(
    well_model::WellModel,
    subsurface_flow_model::GroundwaterFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        flux = check_flux(
            well_model.variables.volumetric_rate[i],
            subsurface_flow_model,
            index,
        )
        well_model.variables.flux[i] = flux
        well_model.variables.flux_av[i] += dt * flux
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

function update_river_storage_stage!(
    gwf_river_model::GwfRiverModel,
    river_flow_model::AbstractRiverFlowModel,
)
    for i in eachindex(gwf_river_model.variables.stage)
        gwf_river_model.variables.stage[i] =
            river_flow_model.variables.h[i] + gwf_river_model.parameters.bottom[i]
        gwf_river_model.variables.storage[i] = river_flow_model.variables.storage[i]
    end
    return nothing
end

update_river_storage_stage!(
    gwf_river_model::Nothing,
    river_flow_model::AbstractRiverFlowModel,
) = nothing

flux!(::Nothing, ::AbstractSubsurfaceFlowModel, ::Vector{Int}, ::Float64) = nothing

get_boundary_index(::RechargeModel, domain::Domain) = domain.land.network.land_indices
get_boundary_index(::GwfRiverModel, domain::Domain) = domain.river.network.land_indices
get_boundary_index(::DrainageModel, domain::Domain) = domain.drain.network.land_indices
get_boundary_index(::Nothing, ::Domain) = Int[]
