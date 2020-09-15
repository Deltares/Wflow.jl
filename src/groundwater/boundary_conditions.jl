struct AquiferBoundaryCondition end


struct River <: AquiferBoundaryCondition {T}
    stage::Vector{T}
    infiltration_conductance::Vector{T}
    exfiltration_conductance::Vector{T}
    bottom::Vector{T}
    index::Vector{Int}
end


function flux!(river::River, aquifer, Q)
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
            
        
struct Drainage <: AquiferBoundaryCondition {T}
    elevation::Vector{T}
    conductance::Vector{T}
    index::Vector{Int}
end


function flux!(drain::Drain, aquifer, Q)
    for (i, index) in enumerate(river.index)
        cond = drain.conductance[i]
        Δϕ = max(0, aquifer.head[index] - drain.elevation[i])
        Q[index] += cond * Δϕ
    end
end


struct HeadBoundary <: AquiferBoundaryCondition {T}
    head::Vector{T}
    conductance::Vector{T}
    index::Vector{Int}
end


function flux!(headboundary::HeadBoundary, Q, aquifer)
    for (i, index) in enumerate(seepage.index)
        cond = headboundary.conductance[i]
        Δϕ = headboundary.head[i] - aquifer.head[index]
        Q[index] += cond * Δϕ
    end
end
