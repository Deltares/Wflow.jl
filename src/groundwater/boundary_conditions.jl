abstract type AquiferBoundaryCondition end


struct River{T} <: AquiferBoundaryCondition
    stage::Vector{T}
    infiltration_conductance::Vector{T}
    exfiltration_conductance::Vector{T}
    bottom::Vector{T}
    index::Vector{Int}
end


function flux!(Q, river::River, aquifer)
    for (i, index) in enumerate(river.index)
        ϕ = aquifer.head[index]
        stage = river.stage[i]
        if stage > ϕ
            cond = river.exfiltration_conductance[i]
            Δϕ = min(stage - bottom, stage - ϕ)
        else
            cond = river.infiltration_conductance[i]
            Δϕ = stage - ϕ
        end
        Q[index] += cond * Δϕ
    end
end
            
        
struct Drainage{T} <: AquiferBoundaryCondition
    elevation::Vector{T}
    conductance::Vector{T}
    index::Vector{Int}
end


function flux!(Q, drain::Drainage, aquifer)
    for (i, index) in enumerate(river.index)
        cond = drain.conductance[i]
        Δϕ = max(0, aquifer.head[index] - drain.elevation[i])
        Q[index] += cond * Δϕ
    end
end


struct HeadBoundary{T} <: AquiferBoundaryCondition
    head::Vector{T}
    conductance::Vector{T}
    index::Vector{Int}
end


function flux!(Q, headboundary::HeadBoundary, aquifer)
    for (i, index) in enumerate(seepage.index)
        cond = headboundary.conductance[i]
        Δϕ = headboundary.head[i] - aquifer.head[index]
        Q[index] += cond * Δϕ
    end
end


struct Recharge{T} <: AquiferBoundaryCondition
    rate::Vector{T}
    index::Vector{Int}
end


function flux!(Q, recharge::Recharge, aquifer)
    for (i, index) in enumerate(recharge.index)
        Q[index] += recharge.rate[i] * aquifer.area[index] 
    end
end


struct Well{T} <: AquiferBoundaryCondition
    volumetric_rate::Vector{T}
    index::Vector{T}
end


function flux!(Q, well::Well, aquifer)
    for (i, index) in enumerate(well.index)
        Q[index] += well.volumetric_rate[i]
    end
end