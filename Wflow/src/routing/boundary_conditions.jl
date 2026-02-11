function check_flux(flux::Float64, gwf::GroundwaterFlow, index::Int)
    # Check if cell is dry
    if gwf.variables.head[index] <= gwf.parameters.bottom[index]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end

function check_flux(flux::Float64, gwf::LateralSSF, index::Int)
    # Check if cell is dry
    if gwf.variables.zi[index] >= gwf.parameters.soilthickness[index]
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
    n::Int
    stage::Vector{Float64} = fill(MISSING_VALUE, n) # [m]
    storage::Vector{Float64} = fill(MISSING_VALUE, n) # [m³]
    flux::Vector{Float64} = fill(MISSING_VALUE, n)  # [m³ d⁻¹]
    flux_av::Vector{Float64} = fill(MISSING_VALUE, n)  # [m³ d⁻¹]
end

@with_kw struct GwfRiver <: AbstractSubsurfaceFlowBC
    parameters::GwfRiverParameters
    variables::GwfRiverVariables
end

function GwfRiver(dataset::NCDataset, config::Config, indices::Vector{CartesianIndex{2}})
    infiltration_conductance = ncread(
        dataset,
        config,
        "river_water__infiltration_conductance";
        optional = false,
        sel = indices,
        type = Float64,
    )
    exfiltration_conductance = ncread(
        dataset,
        config,
        "river_water__exfiltration_conductance";
        optional = false,
        sel = indices,
        type = Float64,
    )
    bottom = ncread(
        dataset,
        config,
        "river_bottom__elevation";
        optional = false,
        sel = indices,
        type = Float64,
    )

    parameters =
        GwfRiverParameters(infiltration_conductance, exfiltration_conductance, bottom)
    n = length(indices)
    variables = GwfRiverVariables(; n)
    river = GwfRiver(parameters, variables)
    return river
end

function flux!(
    river::GwfRiver,
    gwf::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        head = gwf.variables.head[index]
        stage = river.variables.stage[i]
        if stage > head
            max_infiltration_flux = river.variables.storage[i] / dt
            cond = river.parameters.infiltration_conductance[i]
            delta_head = min(stage - river.parameters.bottom[i], stage - head)
            flux = min(cond * delta_head, max_infiltration_flux)
        else
            cond = river.parameters.exfiltration_conductance[i]
            delta_head = stage - head
            flux = check_flux(cond * delta_head, gwf, index)
        end
        river.variables.flux[i] = flux
        gwf.variables.q_net[index] += flux
        river.variables.storage[i] -= dt * flux
        river.variables.flux_av[i] += dt * flux
    end
    return nothing
end

@with_kw struct DrainageParameters
    elevation::Vector{Float64} # [m]
    conductance::Vector{Float64} # [m² d⁻¹]
end

@with_kw struct DrainageVariables
    n::Int
    flux::Vector{Float64} = fill(MISSING_VALUE, n) # [m³ d⁻¹]
    flux_av::Vector{Float64} = fill(MISSING_VALUE, n) # [m³ d⁻¹]
end

@with_kw struct Drainage <: AbstractSubsurfaceFlowBC
    parameters::DrainageParameters
    variables::DrainageVariables
end

function Drainage(dataset::NCDataset, config::Config, indices::Vector{CartesianIndex{2}})
    elevation = ncread(
        dataset,
        config,
        "land_drain__elevation";
        optional = false,
        sel = indices,
        type = Float64,
        fill = MISSING_VALUE,
    )
    conductance = ncread(
        dataset,
        config,
        "land_drain__conductance";
        optional = false,
        sel = indices,
        type = Float64,
        fill = MISSING_VALUE,
    )
    parameters = DrainageParameters(; elevation, conductance)
    n = length(indices)
    variables = DrainageVariables(; n)

    drains = Drainage(parameters, variables)
    return drains
end

function flux!(
    drainage::Drainage,
    gwf::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        cond = drainage.parameters.conductance[i]
        delta_head = min(0, drainage.parameters.elevation[i] - gwf.variables.head[index])
        flux = check_flux(cond * delta_head, gwf, index)
        drainage.variables.flux[i] = flux
        drainage.variables.flux_av[i] += dt * flux
        gwf.variables.q_net[index] += flux
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
    gwf::GroundwaterFlow,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        cond = headboundary.parameters.conductance[i]
        delta_head = headboundary.variables.head[i] - gwf.variables.head[index]
        flux = check_flux(cond * delta_head, gwf, index)
        headboundary.variables.flux[i] = flux
        headboundary.variables.flux_av[i] += dt * flux
        gwf.variables.q_net[index] += flux
    end
    return nothing
end

@with_kw struct RechargeVariables
    n::Int
    rate::Vector{Float64} = fill(MISSING_VALUE, n) # [m d⁻¹]
    flux::Vector{Float64} = zeros(n) # [m³ d⁻¹]
    flux_av::Vector{Float64} = zeros(n) # [m³ d⁻¹]
end

@with_kw struct Recharge <: AbstractSubsurfaceFlowBC
    n::Int
    variables::RechargeVariables = RechargeVariables(; n)
end

function flux!(
    recharge::Recharge,
    gwf::AbstractSubsurfaceFlowModel,
    indices::Vector{Int},
    dt::Float64,
)
    for (i, index) in enumerate(indices)
        flux =
            check_flux(recharge.variables.rate[i] * gwf.parameters.area[index], gwf, index)
        recharge.variables.flux[i] = flux
        recharge.variables.flux_av[i] += dt * flux
        gwf.variables.q_net[index] += flux
    end
    return nothing
end

@with_kw struct WellVariables
    volumetric_rate::Vector{Float64} # [m³ d⁻¹]
    flux::Vector{Float64} # [m³ d⁻¹]
    flux_av::Vector{Float64} # [m³ d⁻¹]
end

@with_kw struct Well <: AbstractSubsurfaceFlowBC
    variables::WellVariables
end

function flux!(well::Well, gwf::GroundwaterFlow, indices::Vector{Int}, dt::Float64)
    for (i, index) in enumerate(indices)
        flux = check_flux(well.variables.volumetric_rate[i], gwf, index)
        well.variables.flux[i] = flux
        well.variables.flux_av[i] += dt * flux
        gwf.variables.q_net[index] += flux
    end
    return nothing
end

function update_river_storage_stage!(river_bc::GwfRiver, river_flow::AbstractRiverFlowModel)
    for i in eachindex(river_bc.variables.stage)
        river_bc.variables.stage[i] =
            river_flow.variables.h[i] + river_bc.parameters.bottom[i]
        river_bc.variables.storage[i] = river_flow.variables.storage[i]
    end
    return nothing
end

update_river_storage_stage!(river_bc::Nothing, river_flow::AbstractRiverFlowModel) = nothing

flux!(::Nothing, ::AbstractSubsurfaceFlowModel, ::Vector{Int}, ::Float64) = nothing

get_boundary_index(::Recharge, domain::Domain) = domain.land.network.land_indices
get_boundary_index(::GwfRiver, domain::Domain) = domain.river.network.land_indices
get_boundary_index(::Drainage, domain::Domain) = domain.drain.network.land_indices
get_boundary_index(::Nothing, ::Domain) = Int[]
