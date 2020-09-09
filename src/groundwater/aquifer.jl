struct Aquifer end


"""
Confined aquifers are overlain by a poorly permeable confining layer (e.g.
clay). No air can get in to fill the pore space so that the aquifer always
remains fully saturated.

For a confined aquifer, water will always flow along the complete height
over the aquifer and transmissivity kH (m d⁻¹) is a constant. Under the
Dupuit-Forchheimer assumption, steady-state groundwater flow in an aquifer is
given by:

Q = kH * dϕ/dx

Specific storage is the amount of water an aquifer releases per unit change in
hydraulic head, per unit volume of aquifer, as the aquifer and the groundwater
itself is compressed. Its value is much smaller than specific yield, between
1E-5 (stiff) and 0.01 (weak).
"""
struct ConfinedAquifer{T,R,L} <: Aquifer
    head::Vector{T}  # hydraulic head [m]
    k::Vector{T}  # horizontal conductivity [m m⁻¹]
    top::Vector{T} # top of groundwater layer [m]
    bottom::Vector{T} # bottom of groundwater layer
    specific_storage::Vector{T} # [m m⁻¹]
    conductance::Vector{T} # Confined aquifer conductance is constant
    area::Vector{T} # area of cell
end


"""
The upper boundary of an unconfined aquifer is the water table (the phreatic
surface).

In an unconfined aquifer, the transmissivity is a function of groundwater head,
resulting in a non-linear equation:

Q = kH(ϕ) * dϕ/dx

Specific yield (or drainable porosity) represents the volumetric fraction the
aquifer will yield when all water drains and the pore volume is filled by air
instead. Specific yield will vary roughly between 0.05 (clay) and 0.45 (peat)
(Johnson, 1967).
"""
struct UnconfinedAquifer{T,R,L} <: Aquifer
    head::Vector{T}  # hydraulic head [m]
    k::Vector{T}  # horizontal conductivity [m d⁻¹]
    top::Vector{T} # top of groundwater layer [m]
    bottom::Vector{T} # bottom of groundwater layer
    specific_yield::Vector{T} # [m m⁻¹]
    # Unconfined aquifer conductance is re-computed
end


"""
The harmonic mean is the exact interblock transmissivity for steady-state
one-dimensional flow with no recharge if the transmissivity is assumed to be
spatially uniform over each finite difference block, changing abruptly at the
block interface. 

Refer to:
    
    Goode, D. J., & Appel, C. A. (1992). Finite-Difference Interblock Transmissivity
    for Unconﬁned Aquifers and for Aquifers having Smoothly Varying Transmissivity.
    Water-resources investigations report, 92, 4124.
"""
function harmonicmean_conductance(k1, k2, D1, D2, l1, l2, width)
    kD1 = k1 * D1
    kD2 = k2 * D2
    if (kD1 * kD2) > 0.0
        return width * kD1 * kD2 / (kD1 * l2 + kD2 * l1)
    else
        return 0.0
    end
end


function saturated_thickness(aquifer::UnconfinedAquifer, head, index)
    min(aquifer.top[index], aquifer.head[index]) - aquifer.bottom[index]
end


function saturated_thickness(aquifer::ConfinedAquifer, head)
    aquifer.top[index] - aquifer.bottom[index]
end

"""
Compute horizontal conductance for a single connection between two cells
(indexed with `i` and `j`).
"""
function horizontal_conductance(
    j::Int,
    i::Int,
    nzi::Int,
    aquifer <: Aquifer,
    M::Mesh,
)
    k1 = aquifer.k[j]
    k2 = aquifer.k[i]
    D1 = saturated_thickness(aquifer, head, j)
    D2 = saturated_thickness(aquifer, head, i)
    length1 = M.length1[nzi]
    length2 = M.length2[nzi]
    width = M.width[nzi]
    return harmonicmean_conductance(k1, k2, D1, D2, length1, length2, width)
end


function darcy_flow!(aquifer::ConfinedAquifer, Q, graph, head, Δt)
    for (nzi, (i, j)) in enumerate(zip(graph))
        ΔH = head[i] - head[j]
        Q[nzi] = -aquifer.conductance[nzi] * ΔH * Δt
    end
end


function darcy_flow!(aquifer::UnconfinedAquifer, Q, M, Δt)
    for j in 1:M.ncell
        # Loop over connections for cell j
        for nzi in connections(M, j)
            # connection from j -> i
            i = M.connection[nzi]
            cond = horizontal_conductance(j, i, nzi, aquifer, graph)
            Q[nzi] = -cond * ΔH * Δt
        end
    end
end


function update(aquifer <: Aquifer, config)
    Δt = Second(config.timestepsecs)
    darcy_flow!(aquifer, Q, connectivity, Δt)
end
