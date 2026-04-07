function check_flux(flux::Float64, subsurface_flow_model::GroundwaterFlowModel, index::Int)
    # Check if cell is dry
    #  [m] <= [m]
    if subsurface_flow_model.variables.head[index] <=
       subsurface_flow_model.parameters.bottom[index]
        # If cell is dry, no negative flux is allowed
        # max([m³ s⁻¹], [m³ s⁻¹])
        return max(0, flux)
    else
        return flux
    end
end

function check_flux(flux::Float64, subsurface_flow_model::LateralSSFModel, index::Int)
    # Check if cell is dry
    # [m] > [m]
    if subsurface_flow_model.variables.zi[index] >=
       subsurface_flow_model.parameters.soilthickness[index]
        # If cell is dry, no negative flux is allowed
        # max([m³ s⁻¹], [m³ s⁻¹])
        return max(0, flux)
    else
        return flux
    end
end

@with_kw struct GwfRiverParameters
    # [m² d⁻¹ => m² s⁻¹]
    infiltration_conductance::Vector{Float64}
    # [m² d⁻¹ => m² s⁻¹]
    exfiltration_conductance::Vector{Float64}
    bottom::Vector{Float64} # [m]
end

@with_kw struct GwfRiverVariables
    n::Int
    # [m]
    stage::Vector{Float64} = fill(MISSING_VALUE, n)
    # [m³]
    storage::Vector{Float64} = fill(MISSING_VALUE, n)
    # [m³ d⁻¹ => m³ s⁻¹]
    flux::Vector{Float64} = fill(MISSING_VALUE, n)
    # [m³ d⁻¹ => m³ s⁻¹]
    flux_av::AverageVector = AverageVector(; n)
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

    parameters =
        GwfRiverParameters(infiltration_conductance, exfiltration_conductance, bottom)
    n = length(indices)
    variables = GwfRiverVariables(; n)
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
        # [m]
        head = subsurface_flow_model.variables.head[index]
        # [m]
        stage = gwf_river_model.variables.stage[i]
        # [m] > [m]
        if stage > head
            # [m³ s⁻¹] = [m³] / [s]
            max_infiltration_flux = gwf_river_model.variables.storage[i] / dt
            # [m² s⁻¹]
            cond = gwf_river_model.parameters.infiltration_conductance[i]
            # [m] = min([m] - [m], [m] - [m])
            delta_head = min(stage - gwf_river_model.parameters.bottom[i], stage - head)
            # [m³ s⁻¹] = min([m² s⁻¹] * [m], [m³ s⁻¹])
            flux = min(cond * delta_head, max_infiltration_flux)
        else
            # [m² s⁻¹]
            cond = gwf_river_model.parameters.exfiltration_conductance[i]
            # [m] = [m] - [m]
            delta_head = stage - head
            # [m³ s⁻¹] = [m² s⁻¹] * [m]
            flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        end
        # [m³ s⁻¹] = [m³ s⁻¹]
        gwf_river_model.variables.flux[i] = flux
        # [m³ s⁻¹] += [m³ s⁻¹]
        subsurface_flow_model.variables.q_net_bnds[index] += flux
        # [m³] -= [s] * [m³ s⁻¹]
        gwf_river_model.variables.storage[i] -= dt * flux
        # [m³] += [s] * [m³ s⁻¹]
        add_to_cumulative!(gwf_river_model.variables.flux_av, i, flux, dt)
    end
    return nothing
end

@with_kw struct DrainageParameters
    # [m]
    elevation::Vector{Float64}
    # [m² d⁻¹ => m² s⁻¹]
    conductance::Vector{Float64}
end

@with_kw struct DrainageVariables
    n::Int
    # [m³ d⁻¹ => m³ s⁻¹]
    flux::Vector{Float64} = fill(MISSING_VALUE, n)
    # [m³ d⁻¹ => m³ s⁻¹]
    flux_av::AverageVector = AverageVector(; n)
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
    elevation = ncread(dataset, config, "land_drain__elevation", Routing; sel=indices)
    conductance = ncread(dataset, config, "land_drain__conductance", Routing; sel=indices)
    parameters = DrainageParameters(; elevation, conductance)
    n = length(indices)
    variables = DrainageVariables(; n)

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
        # [m² s⁻¹]
        cond = drainage_model.parameters.conductance[i]
        # [m] = [m] - [m]
        delta_head = min(
            0,
            drainage_model.parameters.elevation[i] -
            subsurface_flow_model.variables.head[index],
        )
        # [m³ s⁻¹] = [m² s⁻¹] * [m]
        flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        # [m³ s⁻¹] = [m³ s⁻¹]
        drainage_model.variables.flux[i] = flux
        # [m³] += [m³ s⁻¹] * [s]
        add_to_cumulative!(drainage_model.variables.flux_av, i, flux, dt)
        # [m³ s⁻¹] += [m³ s⁻¹]
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@with_kw struct HeadBoundaryParameters
    # [m² d⁻¹] => [m² s⁻¹]
    conductance::Vector{Float64}
end

@with_kw struct HeadBoundaryVariables
    # [m]
    head::Vector{Float64}
    # [m³ d⁻¹ => m³ s⁻¹]
    flux::Vector{Float64}
    # [m³ d⁻¹ => m³ s⁻¹]
    flux_av::AverageVector
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
    for (i, index) in enumerate(indices)
        # [m² s⁻¹]
        cond = headboundary.parameters.conductance[i]
        # [m] = [m] - [m]
        delta_head =
            headboundary.variables.head[i] - subsurface_flow_model.variables.head[index]
        # [m³ s⁻¹] = [m² s⁻¹] * [m]
        flux = check_flux(cond * delta_head, subsurface_flow_model, index)
        # [m³ s⁻¹] = [m³ s⁻¹]
        headboundary.variables.flux[i] = flux
        # [m³] += [m³ s⁻¹] * [s]
        add_to_cumulative!(headboundary.variables.flux_av, i, flux, dt)
        # [m³ s⁻¹] += [m³ s⁻¹]
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@with_kw struct RechargeVariables
    n::Int
    # [m d⁻¹ => m s⁻¹]
    rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # [m³ d⁻¹ => m³ s⁻¹]
    flux::Vector{Float64} = zeros(n)
    # [m³ d⁻¹ => m³ s⁻¹]
    flux_av::AverageVector = AverageVector(; n)
end

@with_kw struct RechargeModel <: AbstractSubsurfaceFlowBC
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
        # [m³ s⁻¹] = [m s⁻¹] * [m²]
        flux = check_flux(
            recharge_model.variables.rate[i] * subsurface_flow_model.parameters.area[index],
            subsurface_flow_model,
            index,
        )
        # [m³ s⁻¹] = [m³ s⁻¹]
        recharge_model.variables.flux[i] = flux
        # [m³] += [m³ s⁻¹] * [s]
        add_to_cumulative!(recharge_model.variables.flux_av, i, flux, dt)
        # [m³ s⁻¹] += [m³ s⁻¹]
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

@with_kw struct WellVariables
    # [m³ d⁻¹ => m³ s⁻¹]
    volumetric_rate::Vector{Float64}
    # [m³ d⁻¹ => m³ s⁻¹]
    flux::Vector{Float64}
    # [m³ d⁻¹ => m³ s⁻¹]
    flux_av::AverageVector
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
    for (i, index) in enumerate(indices)
        # [m³ s⁻¹]
        flux = check_flux(
            well_model.variables.volumetric_rate[i],
            subsurface_flow_model,
            index,
        )
        # [m³ s⁻¹] = [m³ s⁻¹]
        well_model.variables.flux[i] = flux
        add_to_cumulative!(well_model.variables.flux_av, i, flux, dt)
        # [m³ s⁻¹] += [m³ s⁻¹]
        subsurface_flow_model.variables.q_net_bnds[index] += flux
    end
    return nothing
end

function update_river_storage_stage!(
    gwf_river_model::GwfRiverModel,
    river_flow_model::AbstractRiverFlowModel,
)
    for i in eachindex(gwf_river_model.variables.stage)
        # [m] = [m] + [m]
        gwf_river_model.variables.stage[i] =
            river_flow_model.variables.h[i] + gwf_river_model.parameters.bottom[i]
        # [m³] = [m³]
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
