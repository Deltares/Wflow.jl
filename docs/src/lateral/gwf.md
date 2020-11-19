# Groundwater flow

Single layer groundwater flow is defined in the `struct GroundwaterFlow`, which contains the
following fields:

+ aquifer
+ connectivity
+ constanthead
+ boundaries

```julia
Base.@kwdef struct GroundwaterFlow
    aquifer::A where A <: Aquifer
    connectivity::Connectivity
    constanthead::ConstantHead
    boundaries::Vector{B} where B <: AquiferBoundaryCondition
end
```
The fields of `struct GroundwaterFlow` are described in more detail below.

## Aquifer types
Groundwater flow contains the `Abstract type` `Aquifer`, that can be either a confined
(`ConfinedAquifer`) or unconfined (`UnconfinedAquifer`) aquifer. Groundwater flow is solved
forward in time and central in space.

```@docs
Wflow.Aquifer
```

```@docs
Wflow.ConfinedAquifer
```

```@docs
Wflow.UnconfinedAquifer
```
## Connectivity
The connectivity between cells is defined as follows.

```@docs
Wflow.Connectivity
```

## Constant head
Dirichlet boundary conditions can be specified through the field `constanthead`.

## Aquifer boundary conditions

### River
The flux between river and aquifer is calculated using Darcy's law following the approach in
MODFLOW:

```math
    Q_{riv} =  \Bigg\lbrace{C_{i} \,\text{min}(h_{riv} - B_{riv}, h_{riv} - \phi), \,h_{riv} > \phi \atop C_{e} (h_{riv} - \phi) , \,h_{riv} \leq \phi}
```
where ``Q_{riv}`` is the exchange flux from river to aquifer [L``^3`` T``^{-1}``], ``C_i``
[L``^2`` T``^{-1}``] is the river bed infiltration conductance, ``C_e`` [L``^2`` T``^{-1}``]
is the river bed exfiltration conductance, ``B_{riv}`` the bottom of the river bed [L],
``h_{riv}`` is the rive stage [L] and ``\phi`` is the hydraulic head in the river cell [L].

### Drainage

The flux from drains to the aquifer is calculated as follows:

```math
Q_{drain} = C_{drain} \text{min}(0, h_{drain} - \phi)
```

where ``Q_{drain}`` is the exchange flux from drains to aquifer [L``^3`` T``^{-1}``],
``C_{drain}`` [L``^2`` T``^{-1}``] is the drain conductance, ``h_{drain}`` is the drain
elevation [L] and ``\phi`` is the hydraulic head in the cell with drainage [L].

### Recharge
The recharge flux ``Q_{r}`` to the aquifer is calculated as follows:

```math
Q_{r} = R \, A
```
with ``R`` the recharge rate [L T``^{-1}``] and ``A`` the area [L``^2`` ] of the aquifer
cell.

### Head boundary
This boundary is a fixed head with time outside of the model domain, and is generally used
to avoid an unneccessarily extension of the model domain to the location of the fixed
boundary. The flux from the boundary ``Q_{hb}`` [L``^3`` T``^{-1}``] is calculated as
follows:

```math
Q_{hb} = C_{hb} (\phi_{hb} - \phi)
```
with ``C_{hb}`` the conductance of the head boundary [L``^2`` T``^{-1}``], ``\phi_{hb}`` the
head [L] of the head boundary and  ``\phi`` the head of the aquifer cell.

### Well boundary

A volumetric well rate [L``^3`` T``^{-1}``] can be specified as a boundary condition.


!!! note For an unconfined aquifer the boundary fluxes are checked, in case of a dry aquifer
    cell a negative flux is not allowed.

