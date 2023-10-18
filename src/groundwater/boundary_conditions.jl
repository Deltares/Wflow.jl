function check_flux(flux, aquifer::UnconfinedAquifer, index::Int)
    # Check if cell is dry
    if aquifer.head[index] <= aquifer.bottom[index]
        # If cell is dry, no negative flux is allowed
        return max(0, flux)
    else
        return flux
    end
end


# Do nothing for a confined aquifer: aquifer can always provide flux
check_flux(flux, aquifer::ConfinedAquifer, index::Int) = flux


@get_units @exchange @grid_type @grid_location struct River{T} <: AquiferBoundaryCondition
    stage::Vector{T} | "m"
    infiltration_conductance::Vector{T} | "m2 d-1"
    exfiltration_conductance::Vector{T} | "m2 d-1"
    bottom::Vector{T} | "m"
    flux::Vector{T} | "m3 d-1"
    index::Vector{Int} | "-"
end


function flux!(Q, river::River, aquifer)
    for (i, index) in enumerate(river.index)
        ϕ = aquifer.head[index]
        stage = river.stage[i]
        if stage > ϕ
            cond = river.infiltration_conductance[i]
            Δϕ = min(stage - river.bottom[i], stage - ϕ)
        else
            cond = river.exfiltration_conductance[i]
            Δϕ = stage - ϕ
        end
        river.flux[i] = check_flux(cond * Δϕ, aquifer, index)
        Q[index] += river.flux[i]
    end
    return Q
end


@get_units @exchange @grid_type @grid_location struct Drainage{T} <: AquiferBoundaryCondition
    elevation::Vector{T} | "m"
    conductance::Vector{T} | "m2 d-1"
    flux::Vector{T} | "m3 d-1"
    index::Vector{Int} | "-"
end


function flux!(Q, drainage::Drainage, aquifer)
    for (i, index) in enumerate(drainage.index)
        cond = drainage.conductance[i]
        Δϕ = min(0, drainage.elevation[i] - aquifer.head[index])
        drainage.flux[i] = check_flux(cond * Δϕ, aquifer, index)
        Q[index] += drainage.flux[i]
    end
    return Q
end


@get_units struct HeadBoundary{T} <: AquiferBoundaryCondition
    head::Vector{T} | "m"
    conductance::Vector{T} | "m2 d-1"
    flux::Vector{T} | "m3 d-1"
    index::Vector{Int} | "-"
end


function flux!(Q, headboundary::HeadBoundary, aquifer)
    for (i, index) in enumerate(headboundary.index)
        cond = headboundary.conductance[i]
        Δϕ = headboundary.head[i] - aquifer.head[index]
        headboundary.flux[i] = check_flux(cond * Δϕ, aquifer, index)
        Q[index] += headboundary.flux[i]
    end
    return Q
end


@get_units @exchange @grid_type @grid_location struct Recharge{T} <: AquiferBoundaryCondition
    rate::Vector{T} | "m d-1"
    flux::Vector{T} | "m3 d-1"
    index::Vector{Int} | "-"
end


function flux!(Q, recharge::Recharge, aquifer)
    for (i, index) in enumerate(recharge.index)
        recharge.flux[i] =
            check_flux(recharge.rate[i] * aquifer.area[index], aquifer, index)
        Q[index] += recharge.flux[i]
    end
    return Q
end


@get_units struct Well{T} <: AquiferBoundaryCondition
    volumetric_rate::Vector{T} | "m3 d-1"
    flux::Vector{T} | "m3 d-1"
    index::Vector{Int} | "-"
end


function flux!(Q, well::Well, aquifer)
    for (i, index) in enumerate(well.index)
        well.flux[i] = check_flux(well.volumetric_rate[i], aquifer, index)
        Q[index] += well.flux[i]
    end
    return Q
end
