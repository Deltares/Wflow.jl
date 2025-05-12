function check_flux(flux::Float, aquifer::UnconfinedAquifer, index::Int; dt = 1.0)
    # Limit flux, not smaller than negative flux based on aquifer storage
    if flux < 0.0
        max_outflux = -aquifer.variables.storage[index] / dt
        return max(flux, max_outflux)
    else
        return flux
    end
end

# Do nothing for a confined aquifer: aquifer can always provide flux
check_flux(flux::Float, aquifer::ConfinedAquifer, index::Int; dt = 1.0) = flux

@with_kw struct GwfRiverParameters
    infiltration_conductance::Vector{Float} # [m² d⁻¹]
    exfiltration_conductance::Vector{Float} # [m² d⁻¹]
    bottom::Vector{Float} # [m]
end

@with_kw struct GwfRiverVariables
    stage::Vector{Float} # [m]
    storage::Vector{Float} # [m³]
    flux::Vector{Float}  # [m³ d⁻¹]
end

function GwfRiverVariables(n::Int)
    variables = GwfRiverVariables(;
        stage = fill(MISSING_VALUE, n),
        storage = fill(MISSING_VALUE, n),
        flux = fill(MISSING_VALUE, n),
    )
    return variables
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
    lens = lens_input_parameter(config, "river_water__infiltration_conductance")
    infiltration_conductance = ncread(dataset, config, lens; sel = indices, type = Float)

    lens = lens_input_parameter(config, "river_water__exfiltration_conductance")
    exfiltration_conductance = ncread(dataset, config, lens; sel = indices, type = Float)

    lens = lens_input_parameter(config, "river_bottom__elevation")
    bottom = ncread(dataset, config, lens; sel = indices, type = Float)

    parameters =
        GwfRiverParameters(infiltration_conductance, exfiltration_conductance, bottom)
    n = length(indices)
    variables = GwfRiverVariables(n)
    river = GwfRiver(parameters, variables, index)
    return river
end

function flux!(Q::Vector{Float}, river::GwfRiver, aquifer::Aquifer; dt = 1.0)
    for (i, index) in enumerate(river.index)
        head = aquifer.variables.head[index]
        stage = river.variables.stage[i]
        if stage > head
            cond = river.parameters.infiltration_conductance[i]
            delta_head = min(stage - river.parameters.bottom[i], stage - head)
            flux = min(cond * delta_head, river.variables.storage[i] / dt)
        else
            cond = river.parameters.exfiltration_conductance[i]
            delta_head = stage - head
            flux = check_flux(cond * delta_head, aquifer, index; dt)
        end
        river.variables.flux[i] = flux
        Q[index] += flux
    end
    return Q
end

@with_kw struct DrainageParameters
    elevation::Vector{Float} # [m]
    conductance::Vector{Float} # [m² d⁻¹]
end

@with_kw struct DrainageVariables
    flux::Vector{Float} # [m³ d⁻¹]
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
    lens = lens_input_parameter(config, "land_drain__elevation")
    drain_elevation =
        ncread(dataset, config, lens; sel = indices, type = Float, fill = MISSING_VALUE)

    lens = lens_input_parameter(config, "land_drain__conductance")
    drain_conductance =
        ncread(dataset, config, lens; sel = indices, type = Float, fill = MISSING_VALUE)
    elevation = drain_elevation[index]
    conductance = drain_conductance[index]
    parameters = DrainageParameters(; elevation, conductance)
    variables = DrainageVariables(; flux = fill(MISSING_VALUE, length(index)))

    drains = Drainage(parameters, variables, index)
    return drains
end

function flux!(Q::Vector{Float}, drainage::Drainage, aquifer::Aquifer; dt = 1.0)
    for (i, index) in enumerate(drainage.index)
        cond = drainage.parameters.conductance[i]
        delta_head =
            min(0, drainage.parameters.elevation[i] - aquifer.variables.head[index])
        drainage.variables.flux[i] = check_flux(cond * delta_head, aquifer, index; dt)
        Q[index] += drainage.variables.flux[i]
    end
    return Q
end

@with_kw struct HeadBoundaryParameters
    conductance::Vector{Float} # [m² d⁻¹]
end

@with_kw struct HeadBoundaryVariables
    head::Vector{Float} # [m]
    flux::Vector{Float} # [m³ d⁻¹]
end

@with_kw struct HeadBoundary <: AquiferBoundaryCondition
    parameters::HeadBoundaryParameters
    variables::HeadBoundaryVariables
    index::Vector{Int} # [-]
end

function flux!(Q::Vector{Float}, headboundary::HeadBoundary, aquifer::Aquifer; dt = 1.0)
    for (i, index) in enumerate(headboundary.index)
        cond = headboundary.parameters.conductance[i]
        delta_head = headboundary.variables.head[i] - aquifer.variables.head[index]
        headboundary.variables.flux[i] = check_flux(cond * delta_head, aquifer, index; dt)
        Q[index] += headboundary.variables.flux[i]
    end
    return Q
end

@with_kw struct RechargeVariables
    rate::Vector{Float} # [m d⁻¹]
    flux::Vector{Float} # [m³ d⁻¹]
end

@with_kw struct Recharge <: AquiferBoundaryCondition
    variables::RechargeVariables
    index::Vector{Int}  # [-]
end

function Recharge(rate::Vector{Float}, flux::Vector{Float}, index::Vector{Int})
    variables = RechargeVariables(rate, flux)
    recharge = Recharge(variables, index)
    return recharge
end

function flux!(Q::Vector{Float}, recharge::Recharge, aquifer::Aquifer; dt = 1.0)
    for (i, index) in enumerate(recharge.index)
        recharge.variables.flux[i] = check_flux(
            recharge.variables.rate[i] * aquifer.parameters.area[index],
            aquifer,
            index;
            dt,
        )
        Q[index] += recharge.variables.flux[i]
    end
    return Q
end

@with_kw struct WellVariables
    volumetric_rate::Vector{Float} # [m³ d⁻¹]
    flux::Vector{Float} # [m³ d⁻¹]
end

@with_kw struct Well <: AquiferBoundaryCondition
    variables::WellVariables
    index::Vector{Int} # [-]
end

function flux!(Q::Vector{Float}, well::Well, aquifer::Aquifer; dt = 1.0)
    for (i, index) in enumerate(well.index)
        well.variables.flux[i] =
            check_flux(well.variables.volumetric_rate[i], aquifer, index; dt)
        Q[index] += well.variables.flux[i]
    end
    return Q
end
