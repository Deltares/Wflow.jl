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
    if subsurface_flow_model.variables.water_table_depth[index] >=
       subsurface_flow_model.parameters.soil_thickness[index]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end

@with_kw struct GwfRiverParameters
    # [m² s⁻¹]
    infiltration_conductance::Vector{Float64}
    # [m² s⁻¹]
    exfiltration_conductance::Vector{Float64}
    bottom::Vector{Float64} # [m]
end

@with_kw struct GwfRiverVariables
    n_cells::Int
    # [m]
    stage::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # [m³]
    storage::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # [m³ s⁻¹]
    flux::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # [m³]
    flux_cumulative::Vector{Float64} = zeros(n_cells)
    # [m³ s⁻¹]
    flux_average::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

@with_kw struct GwfRiverModel <: AbstractSubsurfaceFlowBC
    parameters::GwfRiverParameters
    variables::GwfRiverVariables
end

function GwfRiverModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    infiltration_conductance = ncread(
        dataset,
        config,
        "river_water__infiltration_conductance",
        Routing;
        sel = indices,
    )
    exfiltration_conductance = ncread(
        dataset,
        config,
        "river_water__exfiltration_conductance",
        Routing;
        sel = indices,
    )
    bottom = ncread(dataset, config, "river_bottom__elevation", Routing; sel = indices)

    parameters =
        GwfRiverParameters(infiltration_conductance, exfiltration_conductance, bottom)
    n_cells = length(indices)
    variables = GwfRiverVariables(; n_cells)
    river_model = GwfRiverModel(parameters, variables)
    return river_model
end

function flux!(
    gwf_river_model::GwfRiverModel,
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (cell_idx, index) in enumerate(indices)
        head = subsurface_flow_model.variables.head[index]
        stage = gwf_river_model.variables.stage[cell_idx]
        if stage > head
            max_infiltration_flux = gwf_river_model.variables.storage[cell_idx] / dt
            cond = gwf_river_model.parameters.infiltration_conductance[cell_idx]
            delta_head =
                min(stage - gwf_river_model.parameters.bottom[cell_idx], stage - head)
            flux = min(cond * delta_head, max_infiltration_flux)
        else
            cond = gwf_river_model.parameters.exfiltration_conductance[cell_idx]
            delta_head = stage - head
            flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        end
        gwf_river_model.variables.flux[cell_idx] = flux
        subsurface_flow_model.variables.q_net_bnds[index] += flux
        gwf_river_model.variables.storage[cell_idx] -= dt * flux
        gwf_river_model.variables.flux_cumulative[cell_idx] += dt * flux
    end
    return nothing
end

@with_kw struct DrainageParameters
    # [m]
    elevation::Vector{Float64}
    # [m² s⁻¹]
    conductance::Vector{Float64}
end

@with_kw struct DrainageVariables
    n_cells::Int
    # [m³ s⁻¹]
    flux::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # [m³ s⁻¹]
    flux_average::Vector{Float64} = zeros(n_cells)
    # [m³]
    flux_cumulative::Vector{Float64} = zeros(n_cells)
end

@with_kw struct DrainageModel <: AbstractSubsurfaceFlowBC
    parameters::DrainageParameters
    variables::DrainageVariables
end

function DrainageModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    elevation = ncread(dataset, config, "land_drain__elevation", Routing; sel = indices)
    conductance = ncread(dataset, config, "land_drain__conductance", Routing; sel = indices)
    parameters = DrainageParameters(; elevation, conductance)
    n_cells = length(indices)
    variables = DrainageVariables(; n_cells)

    drainage_model = DrainageModel(parameters, variables)
    return drainage_model
end

function flux!(
    drainage_model::DrainageModel,
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (cell_idx, index) in enumerate(indices)
        cond = drainage_model.parameters.conductance[cell_idx]
        delta_head = min(
            0,
            drainage_model.parameters.elevation[cell_idx] -
            subsurface_flow_model.variables.head[index],
        )
        flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        drainage_model.variables.flux[cell_idx] = flux
        drainage_model.variables.flux_cumulative[cell_idx] += flux * dt
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@with_kw struct HeadBoundaryParameters
    # [m² s⁻¹]
    conductance::Vector{Float64}
end

@with_kw struct HeadBoundaryVariables
    # [m]
    head::Vector{Float64}
    # [m³ s⁻¹]
    flux::Vector{Float64}
    # [m³]
    flux_cumulative::Vector{Float64}
    # [m³ s⁻¹]
    flux_average::Vector{Float64}
end

@with_kw struct HeadBoundary <: AbstractSubsurfaceFlowBC
    parameters::HeadBoundaryParameters
    variables::HeadBoundaryVariables
end

function flux!(
    headboundary::HeadBoundary,
    subsurface_flow_model::GroundwaterFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (cell_idx, index) in enumerate(indices)
        cond = headboundary.parameters.conductance[cell_idx]
        delta_head =
            headboundary.variables.head[cell_idx] -
            subsurface_flow_model.variables.head[index]
        flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        headboundary.variables.flux[cell_idx] = flux
        headboundary.variables.flux_cumulative[cell_idx] += flux * dt
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@with_kw struct RechargeVariables
    n_cells::Int
    # [m s⁻¹]
    rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # [m³ s⁻¹]
    flux::Vector{Float64} = zeros(n_cells)
    flux_cumulative::Vector{Float64} = zeros(n_cells)
    # [m³ s⁻¹]
    flux_average::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

@with_kw struct RechargeModel <: AbstractSubsurfaceFlowBC
    n_cells::Int
    variables::RechargeVariables = RechargeVariables(; n_cells)
end

function flux!(
    recharge_model::RechargeModel,
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (cell_idx, index) in enumerate(indices)
        flux = check_flux(
            recharge_model.variables.rate[cell_idx] *
            subsurface_flow_model.parameters.area[index],
            subsurface_flow_model,
            index,
        )
        recharge_model.variables.flux[cell_idx] = flux
        recharge_model.variables.flux_cumulative[cell_idx] += flux * dt
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@with_kw struct WellVariables
    # [m³ s⁻¹]
    volumetric_rate::Vector{Float64}
    # [m³ s⁻¹]
    flux::Vector{Float64}
    # [m³]
    flux_cumulative::Vector{Float64}
    # [m³ s⁻¹]
    flux_average::Vector{Float64}
end

@with_kw struct WellModel <: AbstractSubsurfaceFlowBC
    variables::WellVariables
end

function flux!(
    well_model::WellModel,
    subsurface_flow_model::GroundwaterFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (cell_idx, index) in enumerate(indices)
        flux = check_flux(
            well_model.variables.volumetric_rate[cell_idx],
            subsurface_flow_model,
            index,
        )
        well_model.variables.flux[cell_idx] = flux
        well_model.variables.flux_cumulative[cell_idx] += flux * dt
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

function update_river_storage_stage!(
    gwf_river_model::GwfRiverModel,
    river_flow_model::AbstractRiverFlowModel,
)
    for cell_idx in eachindex(gwf_river_model.variables.stage)
        gwf_river_model.variables.stage[cell_idx] =
            river_flow_model.variables.h[cell_idx] +
            gwf_river_model.parameters.bottom[cell_idx]
        gwf_river_model.variables.storage[cell_idx] =
            river_flow_model.variables.storage[cell_idx]
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
