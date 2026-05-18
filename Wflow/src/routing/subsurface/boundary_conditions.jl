function check_flux(
    flux::Float64,
    subsurface_flow_model::GroundwaterFlowModel,
    cell_idx::Int,
)
    # Check if cell is dry
    if subsurface_flow_model.variables.head[cell_idx] <=
       subsurface_flow_model.parameters.bottom[cell_idx]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end

function check_flux(flux::Float64, subsurface_flow_model::LateralSSFModel, cell_idx::Int)
    # Check if cell is dry
    if subsurface_flow_model.variables.zi[cell_idx] >=
       subsurface_flow_model.parameters.soilthickness[cell_idx]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end

@with_kw struct GwfRiverParameters
    infiltration_conductance::Vector{Float64} # [m² d⁻¹]
    exfiltration_conductance::Vector{Float64} # [m² d⁻¹]
    bottom::Vector{Float64} # [m]
end

@with_kw struct GwfRiverVariables
    n_river_cells::Int
    stage::Vector{Float64} = fill(MISSING_VALUE, n_river_cells) # [m]
    storage::Vector{Float64} = fill(MISSING_VALUE, n_river_cells) # [m³]
    flux::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)  # [m³ d⁻¹]
    flux_av::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)  # [m³ d⁻¹]
end

@with_kw struct GwfRiverModel <: AbstractSubsurfaceFlowBC
    parameters::GwfRiverParameters
    variables::GwfRiverVariables
end

function GwfRiverModel(
    dataset::NCDataset,
    config::Config,
    river_indices_2d::Vector{CartesianIndex{2}},
)
    infiltration_conductance = ncread(
        dataset,
        config,
        "river_water__infiltration_conductance",
        Routing;
        sel = river_indices_2d,
    )
    exfiltration_conductance = ncread(
        dataset,
        config,
        "river_water__exfiltration_conductance",
        Routing;
        sel = river_indices_2d,
    )
    bottom =
        ncread(dataset, config, "river_bottom__elevation", Routing; sel = river_indices_2d)

    parameters =
        GwfRiverParameters(infiltration_conductance, exfiltration_conductance, bottom)
    n_river_cells = length(river_indices_2d)
    variables = GwfRiverVariables(; n_river_cells)
    river_model = GwfRiverModel(parameters, variables)
    return river_model
end

function flux!(
    gwf_river_model::GwfRiverModel,
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (river_cell_idx, cell_idx) in enumerate(indices)
        head = subsurface_flow_model.variables.head[cell_idx]
        stage = gwf_river_model.variables.stage[river_cell_idx]
        if stage > head
            max_infiltration_flux = gwf_river_model.variables.storage[river_cell_idx] / dt
            cond = gwf_river_model.parameters.infiltration_conductance[river_cell_idx]
            delta_head =
                min(stage - gwf_river_model.parameters.bottom[river_cell_idx], stage - head)
            flux = min(cond * delta_head, max_infiltration_flux)
        else
            cond = gwf_river_model.parameters.exfiltration_conductance[river_cell_idx]
            delta_head = stage - head
            flux = check_flux(cond * delta_head, subsurface_flow_model, cell_idx)
        end
        gwf_river_model.variables.flux[river_cell_idx] = flux
        subsurface_flow_model.variables.q_net_bnds[cell_idx] += flux
        gwf_river_model.variables.storage[river_cell_idx] -= dt * flux
        gwf_river_model.variables.flux_av[river_cell_idx] += dt * flux
    end
    return nothing
end

@with_kw struct DrainageParameters
    elevation::Vector{Float64} # [m]
    conductance::Vector{Float64} # [m² d⁻¹]
end

@with_kw struct DrainageVariables
    n_cells::Int
    flux::Vector{Float64} = fill(MISSING_VALUE, n_cells) # [m³ d⁻¹]
    flux_av::Vector{Float64} = fill(MISSING_VALUE, n_cells) # [m³ d⁻¹]
end

@with_kw struct DrainageModel <: AbstractSubsurfaceFlowBC
    parameters::DrainageParameters
    variables::DrainageVariables
end

function DrainageModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
)
    elevation =
        ncread(dataset, config, "land_drain__elevation", Routing; sel = land_indices_2d)
    conductance =
        ncread(dataset, config, "land_drain__conductance", Routing; sel = land_indices_2d)
    parameters = DrainageParameters(; elevation, conductance)
    n_cells = length(land_indices_2d)
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
    for (boundary_idx, cell_idx) in enumerate(indices)
        cond = drainage_model.parameters.conductance[boundary_idx]
        delta_head = min(
            0,
            drainage_model.parameters.elevation[boundary_idx] -
            subsurface_flow_model.variables.head[cell_idx],
        )
        flux = check_flux(cond * delta_head, subsurface_flow_model, cell_idx)
        drainage_model.variables.flux[boundary_idx] = flux
        drainage_model.variables.flux_av[boundary_idx] += dt * flux
        subsurface_flow_model.variables.q_net_bnds[cell_idx] += flux
    end
    return nothing
end

@with_kw struct HeadBoundaryParameters
    conductance::Vector{Float64} # [m² d⁻¹]
end

@with_kw struct HeadBoundaryVariables
    head::Vector{Float64} # [m]
    flux::Vector{Float64} # [m³ d⁻¹]
    flux_av::Vector{Float64} # [m³ d⁻¹]
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
    for (boundary_idx, cell_idx) in enumerate(indices)
        cond = headboundary.parameters.conductance[boundary_idx]
        delta_head =
            headboundary.variables.head[boundary_idx] -
            subsurface_flow_model.variables.head[cell_idx]
        flux = check_flux(cond * delta_head, subsurface_flow_model, cell_idx)
        headboundary.variables.flux[boundary_idx] = flux
        headboundary.variables.flux_av[boundary_idx] += dt * flux
        subsurface_flow_model.variables.q_net_bnds[cell_idx] += flux
    end
    return nothing
end

@with_kw struct RechargeVariables
    n_cells::Int
    rate::Vector{Float64} = fill(MISSING_VALUE, n_cells) # [m d⁻¹]
    flux::Vector{Float64} = zeros(n_cells) # [m³ d⁻¹]
    flux_av::Vector{Float64} = zeros(n_cells) # [m³ d⁻¹]
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
    for (boundary_idx, cell_idx) in enumerate(indices)
        flux = check_flux(
            recharge_model.variables.rate[boundary_idx] *
            subsurface_flow_model.parameters.area[cell_idx],
            subsurface_flow_model,
            cell_idx,
        )
        recharge_model.variables.flux[boundary_idx] = flux
        recharge_model.variables.flux_av[boundary_idx] += dt * flux
        subsurface_flow_model.variables.q_net_bnds[cell_idx] += flux
    end
    return nothing
end

@with_kw struct WellVariables
    volumetric_rate::Vector{Float64} # [m³ d⁻¹]
    flux::Vector{Float64} # [m³ d⁻¹]
    flux_av::Vector{Float64} # [m³ d⁻¹]
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
    for (boundary_idx, cell_idx) in enumerate(indices)
        flux = check_flux(
            well_model.variables.volumetric_rate[boundary_idx],
            subsurface_flow_model,
            cell_idx,
        )
        well_model.variables.flux[boundary_idx] = flux
        well_model.variables.flux_av[boundary_idx] += dt * flux
        subsurface_flow_model.variables.q_net_bnds[cell_idx] += flux
    end
    return nothing
end

function update_river_storage_stage!(
    gwf_river_model::GwfRiverModel,
    river_flow_model::AbstractRiverFlowModel,
)
    for river_cell_idx in eachindex(gwf_river_model.variables.stage)
        gwf_river_model.variables.stage[river_cell_idx] =
            river_flow_model.variables.h[river_cell_idx] +
            gwf_river_model.parameters.bottom[river_cell_idx]
        gwf_river_model.variables.storage[river_cell_idx] =
            river_flow_model.variables.storage[river_cell_idx]
    end
    return nothing
end

update_river_storage_stage!(
    gwf_river_model::Nothing,
    river_flow_model::AbstractRiverFlowModel,
) = nothing

flux!(::Nothing, ::AbstractSubsurfaceFlowModel, ::Vector{Int}, ::Float64) = nothing

get_boundary_index(::RechargeModel, domain::Domain) = domain.land.network.cell_indices
get_boundary_index(::GwfRiverModel, domain::Domain) =
    domain.river.network.cell_indices_containing_river
get_boundary_index(::DrainageModel, domain::Domain) =
    domain.drain.network.cell_indices_containing_drainage
get_boundary_index(::Nothing, ::Domain) = Int[]
