function check_flux(flux, aquifer::UnconfinedAquifer, index::Int)
    # Check if cell is dry
    if aquifer.variables.head[index] <= aquifer.parameters.bottom[index]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end

# Do nothing for a confined aquifer: aquifer can always provide flux
check_flux(flux, aquifer::ConfinedAquifer, index::Int) = flux

@get_units @grid_loc @with_kw struct RiverParameters{T}
    infiltration_conductance::Vector{T} | "m2 d-1"
    exfiltration_conductance::Vector{T} | "m2 d-1"
    bottom::Vector{T} | "m"
end

@get_units @grid_loc @with_kw struct RiverVariables{T}
    stage::Vector{T} | "m"
    flux::Vector{T} | "m3 d-1"
end

function RiverVariables(n)
    variables = RiverVariables{Float}(; stage = fill(mv, n), flux = fill(mv, n))
    return variables
end

@get_units @grid_loc @with_kw struct River{T} <: AquiferBoundaryCondition
    parameters::RiverParameters{T}
    variables::RiverVariables{T}
    index::Vector{Int} | "-"
end

function River(dataset, config, indices, index)
    infiltration_conductance = ncread(
        dataset,
        config,
        "lateral.subsurface.infiltration_conductance";
        sel = indices,
        type = Float,
    )
    exfiltration_conductance = ncread(
        dataset,
        config,
        "lateral.subsurface.exfiltration_conductance";
        sel = indices,
        type = Float,
    )
    bottom = ncread(
        dataset,
        config,
        "lateral.subsurface.river_bottom";
        sel = indices,
        type = Float,
    )

    parameters =
        RiverParameters{Float}(infiltration_conductance, exfiltration_conductance, bottom)
    n = length(indices)
    variables = RiverVariables(n)
    river = River(parameters, variables, index)
    return river
end

function flux!(Q, river::River, aquifer)
    for (i, index) in enumerate(river.index)
        head = aquifer.variables.head[index]
        stage = river.variables.stage[i]
        if stage > head
            cond = river.parameters.infiltration_conductance[i]
            delta_head = min(stage - river.parameters.bottom[i], stage - head)
        else
            cond = river.parameters.exfiltration_conductance[i]
            delta_head = stage - head
        end
        river.variables.flux[i] = check_flux(cond * delta_head, aquifer, index)
        Q[index] += river.variables.flux[i]
    end
    return Q
end

@get_units @grid_loc @with_kw struct DrainageParameters{T}
    elevation::Vector{T} | "m"
    conductance::Vector{T} | "m2 d-1"
end

@get_units @grid_loc @with_kw struct DrainageVariables{T}
    flux::Vector{T} | "m3 d-1"
end

@get_units @grid_loc @with_kw struct Drainage{T} <: AquiferBoundaryCondition
    parameters::DrainageParameters{T}
    variables::DrainageVariables{T}
    index::Vector{Int} | "-"
end

function Drainage(dataset, config, indices, index)
    drain_elevation = ncread(
        dataset,
        config,
        "lateral.subsurface.drain_elevation";
        sel = indices,
        type = Float,
        fill = mv,
    )
    drain_conductance = ncread(
        dataset,
        config,
        "lateral.subsurface.drain_conductance";
        sel = indices,
        type = Float,
        fill = mv,
    )
    elevation = drain_elevation[index]
    conductance = drain_conductance[index]
    parameters = DrainageParameters{Float}(; elevation, conductance)
    variables = DrainageVariables{Float}(; flux = fill(mv, length(index)))

    drains = Drainage{Float}(parameters, variables, index)
    return drains
end

function flux!(Q, drainage::Drainage, aquifer)
    for (i, index) in enumerate(drainage.index)
        cond = drainage.parameters.conductance[i]
        delta_head =
            min(0, drainage.parameters.elevation[i] - aquifer.variables.head[index])
        drainage.variables.flux[i] = check_flux(cond * delta_head, aquifer, index)
        Q[index] += drainage.variables.flux[i]
    end
    return Q
end

@get_units @grid_loc @with_kw struct HeadBoundaryParameters{T}
    conductance::Vector{T} | "m2 d-1"
end

@get_units @grid_loc @with_kw struct HeadBoundaryVariables{T}
    head::Vector{T} | "m"
    flux::Vector{T} | "m3 d-1"
end

@get_units @grid_loc @with_kw struct HeadBoundary{T} <: AquiferBoundaryCondition
    parameters::HeadBoundaryParameters{T}
    variables::HeadBoundaryVariables{T}
    index::Vector{Int} | "-"
end

function flux!(Q, headboundary::HeadBoundary, aquifer)
    for (i, index) in enumerate(headboundary.index)
        cond = headboundary.parameters.conductance[i]
        delta_head = headboundary.variables.head[i] - aquifer.variables.head[index]
        headboundary.variables.flux[i] = check_flux(cond * delta_head, aquifer, index)
        Q[index] += headboundary.variables.flux[i]
    end
    return Q
end

@get_units @grid_loc @with_kw struct RechargeVariables{T}
    rate::Vector{T} | "m d-1"
    flux::Vector{T} | "m3 d-1"
end

@get_units @grid_loc @with_kw struct Recharge{T} <: AquiferBoundaryCondition
    variables::RechargeVariables{T}
    index::Vector{Int} | "-"
end

function Recharge(rate, flux, index)
    variables = RechargeVariables{Float}(rate, flux)
    recharge = Recharge{Float}(variables, index)
    return recharge
end

function flux!(Q, recharge::Recharge, aquifer)
    for (i, index) in enumerate(recharge.index)
        recharge.variables.flux[i] = check_flux(
            recharge.variables.rate[i] * aquifer.parameters.area[index],
            aquifer,
            index,
        )
        Q[index] += recharge.variables.flux[i]
    end
    return Q
end

@get_units @grid_loc @with_kw struct WellVariables{T}
    volumetric_rate::Vector{T} | "m3 d-1"
    flux::Vector{T} | "m3 d-1"
end

@get_units @grid_loc @with_kw struct Well{T} <: AquiferBoundaryCondition
    variables::WellVariables{T}
    index::Vector{Int} | "-"
end

function flux!(Q, well::Well, aquifer)
    for (i, index) in enumerate(well.index)
        well.variables.flux[i] =
            check_flux(well.variables.volumetric_rate[i], aquifer, index)
        Q[index] += well.variables.flux[i]
    end
    return Q
end
