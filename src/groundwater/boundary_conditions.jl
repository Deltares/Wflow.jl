struct AquiferBoundaryCondition end


struct River <: AquiferBoundaryCondition {T,R,L}
    stage::Vector{T}
    infiltration_conductance::Vector{T}
    exfiltration_conductance::Vector{T}
    bottom::Vector{T}
    index::Vector{Int}
end


function flux!(river::River, Q, aquifer)
    for (i, index) in enumerate(river.index)
        ϕ = aquifer.head[index]
        stage = river.stage[i]
        if stage > ϕ
            cond = river.exfiltration_conductance[i]
            dϕ = min(stage - bottom, stage - ϕ)
        else
            cond = river.infiltration_conductance[i]
            dϕ = stage - ϕ
        end
        Q[index] += cond * dϕ
    end
end
            
        
struct Drainage <: AquiferBoundaryCondition {T,R,L}
    elevation::Vector{T}
    conductance::Vector{T}
    index::Vector{Int}
end


function flux!(drain::Drain, Q, aquifer)
    for (i, index) in enumerate(river.index)
        cond = drain.conductance[i]
        dϕ = max(0, aquifer.head[index] - drain.elevation[i])
        Q[index] += cond * dϕ
    end
end


struct HeadBoundary <: AquiferBoundaryCondition {T,R,L}
    head::Vector{T}
    conductance::Vector{T}
    index::Vector{Int}
end


function flux!(headboundary::HeadBoundary, Q, aquifer)
    for (i, index) in enumerate(seepage.index)
        cond = headboundary.conductance[i]
        dϕ = headboundary.head[i] - aquifer.head[index]
        Q[index] += cond * dϕ
    end
end
