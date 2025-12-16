function check_flux(flux::Float64, aquifer::UnconfinedAquifer, index::Int)
    # Check if cell is dry
    if aquifer.variables.head[index] <= aquifer.parameters.bottom[index]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end

# Do nothing for a confined aquifer: aquifer can always provide flux
check_flux(flux::Float64, aquifer::ConfinedAquifer, index::Int) = flux

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

@with_kw struct GwfRiver <: AquiferBoundaryCondition
    parameters::GwfRiverParameters
    variables::GwfRiverVariables
    index::Vector{Int} # [-]
end

function GwfRiver(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    index::Vector{Int},
)
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
    river = GwfRiver(parameters, variables, index)
    return river
end

function flux!(river::GwfRiver, aquifer::Aquifer, dt::Float64)
    for (i, index) in enumerate(river.index)
        head = aquifer.variables.head[index]
        stage = river.variables.stage[i]
        if stage > head
            max_infiltration_flux = river.variables.storage[i] / dt
            cond = river.parameters.infiltration_conductance[i]
            delta_head = min(stage - river.parameters.bottom[i], stage - head)
            flux = min(cond * delta_head, max_infiltration_flux)
        else
            cond = river.parameters.exfiltration_conductance[i]
            delta_head = stage - head
            flux = check_flux(cond * delta_head, aquifer, index)
        end
        river.variables.flux[i] = flux
        aquifer.variables.q_net[index] += flux
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

@with_kw struct Drainage <: AquiferBoundaryCondition
    parameters::DrainageParameters
    variables::DrainageVariables
    index::Vector{Int} # [-]
end

function Drainage(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    index::Vector{Int},
)
    drain_elevation = ncread(
        dataset,
        config,
        "land_drain__elevation";
        optional = false,
        sel = indices,
        type = Float64,
        fill = MISSING_VALUE,
    )
    drain_conductance = ncread(
        dataset,
        config,
        "land_drain__conductance";
        optional = false,
        sel = indices,
        type = Float64,
        fill = MISSING_VALUE,
    )
    elevation = drain_elevation[index]
    conductance = drain_conductance[index]
    parameters = DrainageParameters(; elevation, conductance)
    n = length(index)
    variables = DrainageVariables(; n)

    drains = Drainage(parameters, variables, index)
    return drains
end

function flux!(drainage::Drainage, aquifer::Aquifer, dt::Float64)
    for (i, index) in enumerate(drainage.index)
        cond = drainage.parameters.conductance[i]
        delta_head =
            min(0, drainage.parameters.elevation[i] - aquifer.variables.head[index])
        flux = check_flux(cond * delta_head, aquifer, index)
        drainage.variables.flux[i] = flux
        drainage.variables.flux_av[i] += dt * flux
        aquifer.variables.q_net[index] += flux
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

@with_kw struct HeadBoundary <: AquiferBoundaryCondition
    parameters::HeadBoundaryParameters
    variables::HeadBoundaryVariables
    index::Vector{Int} # [-]
end

function flux!(headboundary::HeadBoundary, aquifer::Aquifer, dt::Float64)
    for (i, index) in enumerate(headboundary.index)
        cond = headboundary.parameters.conductance[i]
        delta_head = headboundary.variables.head[i] - aquifer.variables.head[index]
        flux = check_flux(cond * delta_head, aquifer, index)
        headboundary.variables.flux[i] = flux
        headboundary.variables.flux_av[i] += dt * flux
        aquifer.variables.q_net[index] += flux
    end
    return nothing
end

@with_kw struct RechargeVariables
    n::Int
    rate::Vector{Float64} = fill(MISSING_VALUE, n) # [m d⁻¹]
    flux::Vector{Float64} = zeros(n) # [m³ d⁻¹]
    flux_av::Vector{Float64} = zeros(n) # [m³ d⁻¹]
end

@with_kw struct Recharge <: AquiferBoundaryCondition
    n::Int
    variables::RechargeVariables = RechargeVariables(; n)
    index::Vector{Int} = collect(1:n)  # [-]
end

function flux!(recharge::Recharge, aquifer::Aquifer, dt::Float64)
    for (i, index) in enumerate(recharge.index)
        flux = check_flux(
            recharge.variables.rate[i] * aquifer.parameters.area[index],
            aquifer,
            index,
        )
        recharge.variables.flux[i] = flux
        recharge.variables.flux_av[i] += dt * flux
        aquifer.variables.q_net[index] += flux
    end
    return nothing
end

@with_kw struct WellVariables
    volumetric_rate::Vector{Float64} # [m³ d⁻¹]
    flux::Vector{Float64} # [m³ d⁻¹]
    flux_av::Vector{Float64} # [m³ d⁻¹]
end

@with_kw struct Well <: AquiferBoundaryCondition
    variables::WellVariables
    index::Vector{Int} # [-]
end

function flux!(well::Well, aquifer::Aquifer, dt::Float64)
    for (i, index) in enumerate(well.index)
        flux = check_flux(well.variables.volumetric_rate[i], aquifer, index)
        well.variables.flux[i] = flux
        well.variables.flux_av[i] += dt * flux
        aquifer.variables.q_net[index] += flux
    end
    return nothing
end

flux!(::Nothing, ::Aquifer, ::Float64) = nothing
