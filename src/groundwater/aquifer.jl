"""
The vertically averaged governing equation of an unconfined, inhomogeneous and
isotropic aquifer in one dimension can be written as:

S * ∂ϕ / ∂t = ∂ / ∂x * (kH * ∂ϕ / ∂x) + Q

with:
* S: storativity (or storage coefficient)
* ϕ: hydraulic head (groundwater level)
* t: time
* k: conductivity
* H: H the (saturated) aquifer height: groundwater level - aquifer bottom
  elevation
* η: elevation of aquifer bottom
* Q: fluxes from boundary conditions (e.g. recharge or abstraction)

The simplest finite difference formulation is forward in time, central in space,
and can be written as:

Sᵢ * (ϕᵢᵗ⁺¹ - ϕᵢᵗ) / Δt = -Cᵢ₋₁ * (ϕᵢ₋₁ - ϕᵢ) -Cᵢ * (ϕᵢ₊₁ - ϕᵢ) + Qᵢ

with:
* ᵢ as cell index
* ᵗ as time index
* Δt as step size
* Cᵢ₋₁ as the intercell conductance between cell i-1 and i
* Cᵢ as the intercell conductance between cell i and i+1

Conductance is defined as:

C = kH * w / l

with:
* w the width of the cell to cell connection
* l the length of the cell to cell connection

k and H may both vary in space; intercell conductance is therefore an average
using the properties of two cells. See the documentation below.

There is only one unknown, ϕᵢᵗ⁺¹. Reshuffling terms: 

ϕᵢᵗ⁺¹ = ϕᵢᵗ + (Cᵢ₋₁ * (ϕᵢ - ϕᵢ₋₁) + Cᵢ * (ϕᵢ₊₁ - ϕᵢ) + Qᵢ) * Δt / Sᵢ

This can be generalized to two dimensions, for both regular and irregular cell
connectivity.

See this paper for more details:
     Chu, W. S., & Willis, R. (1984). An explicit finite difference model for
     unconfined aquifers. Groundwater, 22(6), 728-734. 

Boundary conditions can be classified into three categories:
* specified head (Dirichlet)
* specified flux (Neumann)
* head-dependent flux (Robin)

Neumann and Robin conditions are implemented by adding to or subtracting from a
net (lumped) cell flux. Dirichlet conditions are special cased, since they
cannot (easily) be implemented via the flux, but the head is set directly
instead.
"""
abstract type Aquifer end


"""
Confined aquifers are overlain by a poorly permeable confining layer (e.g.
clay). No air can get in to fill the pore space so that the aquifer always
remains fully saturated. For a confined aquifer, water will always flow along
the complete height over the aquifer and transmissivity kH (m d⁻¹) is a
constant.

Specific storage is the amount of water an aquifer releases per unit change in
hydraulic head, per unit volume of aquifer, as the aquifer and the groundwater
itself is compressed. Its value is much smaller than specific yield, between
1E-5 (stiff) and 0.01 (weak).

NOTA BENE: **specific** storage is per m of aquifer (conf. specific weight).
**Storativity** or (**storage coefficient**) is for the entire aquifer (conf.
transmissivity).
"""
struct ConfinedAquifer{T} <: Aquifer
    head::Vector{T}  # hydraulic head [m]
    k::Vector{T}  # horizontal conductivity [m m⁻¹]
    top::Vector{T} # top of groundwater layer [m]
    bottom::Vector{T} # bottom of groundwater layer
    area::Vector{T} # area of cell
    specific_storage::Vector{T} # [m m⁻¹ m⁻¹]
    storativity::Vector{T} # [m m⁻¹]
    conductance::Vector{T} # Confined aquifer conductance is constant
end


"""
The upper boundary of an unconfined aquifer is the water table (the phreatic
surface).

Specific yield (or drainable porosity) represents the volumetric fraction the
aquifer will yield when all water drains and the pore volume is filled by air
instead. Specific yield will vary roughly between 0.05 (clay) and 0.45 (peat)
(Johnson, 1967).
"""
struct UnconfinedAquifer{T} <: Aquifer
    head::Vector{T}  # hydraulic head [m]
    k::Vector{T}  # horizontal conductivity [m d⁻¹]
    top::Vector{T} # top of groundwater layer [m]
    bottom::Vector{T} # bottom of groundwater layer
    area::Vector{T}
    specific_yield::Vector{T} # [m m⁻¹]
    # Unconfined aquifer conductance is re-computed
end


storativity(A::UnconfinedAquifer) = A.specific_yield
storativity(A::ConfinedAquifer) = A.storativity


"""
    harmonicmean_conductance(k1, k2, H1, H2, l1, l2, width)

The harmonic mean is the exact interblock transmissivity for steady-state
one-dimensional flow with no recharge if the transmissivity is assumed to be
spatially uniform over each finite difference block, changing abruptly at the
block interface. 

Refer to:
    
    Goode, D. J., & Appel, C. A. (1992). Finite-Difference Interblock Transmissivity
    for Unconﬁned Aquifers and for Aquifers having Smoothly Varying Transmissivity.
    Water-resources investigations report, 92, 4124.
"""
function harmonicmean_conductance(k1, k2, H1, H2, l1, l2, width)
    kH1 = k1 * H1
    kH2 = k2 * H2
    if (kH1 * kH2) > 0.0
        return width * kH1 * kH2 / (kH1 * l2 + kH2 * l1)
    else
        return 0.0
    end
end


function saturated_thickness(aquifer::UnconfinedAquifer, index::Int)
    min(aquifer.top[index], aquifer.head[index]) - aquifer.bottom[index]
end


function saturated_thickness(aquifer::ConfinedAquifer, index::Int)
    aquifer.top[index] - aquifer.bottom[index]
end


"""
    horizontal_conductance(i, j, nzi, aquifer, C)

Compute horizontal conductance for a single connection between two cells
(indexed with `i` and `j`). Geometry characteristics are taken from the
connectivity struct `C`, using the non-zero index (nzi) of its CSC data
structure.
"""
function horizontal_conductance(
    i::Int,
    j::Int,
    nzi::Int,
    aquifer::A,
    connectivity::Connectivity,
) where A <: Aquifer
    k1 = aquifer.k[i]
    k2 = aquifer.k[j]
    H1 = saturated_thickness(aquifer, i)
    H2 = saturated_thickness(aquifer, j)
    length1 = connectivity.length1[nzi]
    length2 = connectivity.length2[nzi]
    width = connectivity.width[nzi]
    return harmonicmean_conductance(k1, k2, H1, H2, length1, length2, width)
end

"""
Conductance for a confined aquifer is constant, and only has to be set once.
"""
function initialize_conductance!(aquifer::ConfinedAquifer, connectivity::Connectivity)
    for i in 1:connectivity.ncell
        # Loop over connections for cell j
        for nzi in connections(connectivity, i)
            j = connectivity.rowval[nzi]
            aquifer.conductance[nzi] = horizontal_conductance(i, j, nzi, aquifer, connectivity)
        end
    end
end


function conductance(aquifer::ConfinedAquifer, connectivity, i, j, nzi)
    return aquifer.conductance[nzi]
end


function conductance(aquifer::UnconfinedAquifer, connectivity, i, j, nzi)
    return horizontal_conductance(i, j, nzi, aquifer, connectivity)
end


function flux!(Q, aquifer, connectivity)
    for i in 1:connectivity.ncell
        # Loop over connections for cell j
        for nzi in connections(connectivity, i)
            # connection from i -> j
            j = connectivity.rowval[nzi]
            Δϕ = aquifer.head[i] - aquifer.head[j]
            cond = conductance(aquifer, connectivity, i, j, nzi)
Q[i] -= cond * Δϕ
        end
    end
end


struct ConstantHead{T}
    head::Vector{T}
    index::Vector{Int}
end


"""
Compute a stable timestep size given the forward-in-time, central in space scheme.
The following criterion can be found in Chu & Willis (1984)

Δt * k * H / (Δx * Δy * S) <= 1/4
"""
function stable_timestep(aquifer)
    Δtₘᵢₙ = Inf
    for i in eachindex(aquifer.head)
        Δt = aquifer.area[i] * storativity(aquifer)[i] / (aquifer.k[i] * saturated_thickness(aquifer, i))
        Δtₘᵢₙ = Δt < Δtₘᵢₙ ? Δt : Δtₘᵢₙ
    end
    return 0.25 * Δtₘᵢₙ
end


function update(gwf, Q, Δt)
    Q .= 0.0  # TODO: Probably remove this when linking with other components
    flux!(Q, gwf.aquifer, gwf.connectivity)
    for boundary in gwf.boundaries
        flux!(Q, boundary, gwf.aquifer)
    end
    gwf.aquifer.head .+= (Q ./ gwf.aquifer.area .* Δt ./ storativity(gwf.aquifer))
    # Set constant head (dirichlet) boundaries
    gwf.aquifer.head[gwf.constanthead.index] .= gwf.constanthead.head
    # Make sure no heads ends up below the aquifer bottom
    gwf.aquifer.head .= max.(gwf.aquifer.head, gwf.aquifer.bottom)
end


Base.@kwdef struct GroundwaterFlow
    aquifer::A where A <: Aquifer
    connectivity::Connectivity
    constanthead::ConstantHead
    boundaries::Vector{B} where B <: AquiferBoundaryCondition
end 