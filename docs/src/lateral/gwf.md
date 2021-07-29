# [Groundwater flow](@id lateral_gwf)

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
(`ConfinedAquifer`) or unconfined (`UnconfinedAquifer`) aquifer. Confined aquifers are
overlain by a poorly permeable confining layer (e.g. clay). No air can get in to fill the
pore space so that the aquifer always remains fully saturated. For a confined aquifer, water
will always flow along the complete height ``H`` [m] over the aquifer and transmissivity
``kH`` [m``^2`` d``^{-1}``] is a constant (``k`` [m d``^{-1}``] is the horizontal hydraulic
conductivity). Specific storage is the amount of water an aquifer releases per unit change in
hydraulic head, per unit volume of aquifer, as the aquifer and the groundwater itself is
compressed. Its value is much smaller than specific yield, between 1e-5 (stiff) and 0.01
(weak). The upper boundary of an unconfined aquifer is the water table (the phreatic
surface). Specific yield (or drainable porosity) represents the volumetric fraction the
aquifer will yield when all water drains and the pore volume is filled by air instead.
Specific yield will vary roughly between 0.05 (clay) and 0.45 (peat) (Johnson, 1967).

The parameters of structs `ConfinedAquifer` and `UnconfinedAquifer` are presented in Tables
of the [Confined aquifer](@ref) and [Unconfined aquifer](@ref) sections of Model parameters.
Parameters that can be set directly from the static input data (netCDF) are marked in these
Tables.

Groundwater flow is solved forward in time and central in space. The vertically averaged
governing equation for an inhomogeneous and isotropic aquifer in one dimension can be
written as:

```math
    S \frac{\phi}{\delta t} =  \frac{\delta}{\delta x} (kH \frac{\phi}{\delta x}) + Q
```

where ``S`` [m m``^{-1}``] is storativity (or specific yield), ``\phi`` [m] is hydraulic
head, ``t`` is time, ``k`` [m t``^{-1}``] is horizontal hydraulic conductivity, ``H`` [m] is
the (saturated) aquifer height: groundwater level - aquifer bottom elevation and ``Q`` [m
t``^{-1}``] represents fluxes from boundary conditions (e.g. recharge or abstraction), see
also [Aquifer boundary conditions](@ref).

The simplest finite difference formulation is forward in time, central in space, and can be
written as:

```math
    S_i  \frac{(\phi_{i}^{t+1} - \phi_i^{t})}{\Delta t} = -C_{i-1}  (\phi_{i-1} - \phi_i) - C_i  (\phi_{i+1} - \phi_i) + Q_ᵢ
```

where ``_i`` is the cell index, ``^t`` is time, ``\Delta t`` is the step size, ``C_{i-1}``
is the the intercell conductance between cell ``i-1`` and ``i`` and ``C_i`` is the intercell
conductance between cell ``i`` and ``i+1``. The connection data between cells is stored as
part of the `Connectivity` struct, see also [Connectivity](@ref) for more information.

Conductance ``C`` is defined as:

```math 
    C = \frac{kH w}{l}
```

where ``w`` [m] is the width of the cell to cell connection, and ``l`` [m] is the length of
the cell to cell connection. ``k`` and ``H`` may both vary in space; intercell conductance
is therefore an average using the properties of two cells.  For the calculation of the
intercell conductance ``C`` the harmonic mean is used (see also Goode and Appel, 1992), here
between cell index ``i`` and cell index ``i+1``, in the ``x`` direction:

```math
    C_i = w  \frac{(k_iH_i\cdot k_{i+1}H_{i+1})}{(k_iH_i \cdot l_{i+1} + k_{i+1}H_{i+1} \cdot l_i)}
```

where ``H`` [m] is the aquifer top - aquifer bottom, and ``k``, ``l_i`` is the length in
cell ``i`` (``0.5 \Delta x_i``),  ``l_{i+1}`` is the length in cell ``i+1`` (``0.5 \Delta
x_{i+1}``) and ``w`` as previously defined. For an unconfined aquifer the intercell
conductance is scaled by using the "upstream saturated fraction" as the MODFLOW
documentation calls it. In this approach, the saturated thickness of a cell-to-cell is
approximated using the cell with the highest head. This results in a consistent
overestimation of the saturated thickness, but it avoids complexities related with cell
drying and rewetting, such as having to define a "wetting threshold" or a "wetting factor".
See also the documentation for MODFLOW-NWT (Niswonger et al., 2011) or MODFLOW6 (Langevin et
al., 2017) for more background information. For more background on drying and rewetting, see
for example McDonald et al. (1991).

For the finite difference formulation, there is only one unknown, ``\phi_i^{t+1}``.
Reshuffling terms:

```math
\phi_i^{t+1} = \phi_i^t + (C_{i-1}  (\phi_i - \phi_{i-1}) + C_i  (\phi_{i+1} - \phi_i) + Q_i) \frac{Δt}{S_i}
```

This can be generalized to two dimensions, for both regular and irregular cell connectivity.
Finally, a stable time step size can be computed given the forward-in-time, central in space
scheme, based on the following criterion from Chu and Willis (1984):

```math
\frac{\Delta t k H}{(\Delta x  \Delta y S)}  \le \frac{1}{4}
```
where ``\Delta t`` [d] is the stable time step size, ``\Delta x`` [m] is the cell length in
the ``x`` direction and ``\Delta y`` [m] is the cell length in the ``y`` direction, ``k`` is
the horizontal hydraulic conductivity [m``^2`` d``^{-1}``] and ``H`` [m] is the saturated
thickness of the aquifer. For each cell ``\frac{(\Delta x  \Delta y S)}{k H}`` is
calculated, the minimum of these values is determined, and multiplied by ``\frac{1}{4}``, to
get the stable time step size.

For more details about the finite difference formulation and the stable time step size
criterion we refer to the paper of Chu and Willis (1984).

Boundary conditions can be classified into three categories:

+ specified head (Dirichlet)
+ specified flux (Neumann)
+ head-dependent flux (Robin)

Neumann and Robin conditions are implemented by adding to or subtracting from a net (lumped)
cell flux. Dirichlet conditions are special cased, since they cannot (easily) be implemented
via the flux, but the head is set directly instead.

## Connectivity
The connectivity between cells is defined as follows.

```@docs
Wflow.Connectivity
```

## Constant head
Dirichlet boundary conditions can be specified through the field `constanthead` (type
`ConstantHead`)  of the `GroundwaterFlow` struct.

```julia
@get_units struct ConstantHead{T}
    head::Vector{T} | "m"
    index::Vector{Int} | "-"
end
```

For the model `SBM + Groundwater flow` this boundary condition is optional, and if used
should be specified in the TOML file as follows (see also
[sbm\_gwf\_config.toml](https://github.com/Deltares/Wflow.jl/blob/master/test/sbm_gwf_config.toml)):

```toml
[model]
constanthead = true
```

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
``h_{riv}`` is the river stage [L] and ``\phi`` is the hydraulic head in the river cell [L].

The Table in the Groundwater flow [river boundary condition](@ref gwf_river_params) section
of the Model parameters provides the parameters of the struct `River`. Parameters that can
be set directly from the static input data (netCDF) are marked in this Table.

The exchange flux (river to aquifer) ``Q_{riv}`` is an output variable (field `flux` of the
`River` struct), and is used to update the total flux in a river cell. For the model `SBM +
Groundwater flow`, the water level `h` [m] of the river kinematic wave in combination with
the river `bottom` is used to update the `stage` field of the `River` struct each time step.

### Drainage
The flux from drains to the aquifer is calculated as follows:

```math
Q_{drain} = C_{drain} \text{min}(0, h_{drain} - \phi)
```

where ``Q_{drain}`` is the exchange flux from drains to aquifer [L``^3`` T``^{-1}``],
``C_{drain}`` [L``^2`` T``^{-1}``] is the drain conductance, ``h_{drain}`` is the drain
elevation [L] and ``\phi`` is the hydraulic head in the cell with drainage [L].

The Table in the Groundwater flow [drainage boundary condition](@ref gwf_drainage_params)
section of the Model parameters provides the parameters of the struct `Drainage`. Parameters
that can be set directly from the static input data (netCDF) are marked in this Table.

The exchange flux (drains to aquifer) ``Q_{drain}`` is an output variable (field `flux` of
struct `Drainage`), and is used to update the total flux in a cell with drains. For the
model `SBM + Groundwater flow` this boundary condition is optional, and if used should be
specified in the TOML file as follows (see also
[sbm\_gwf\_config.toml](https://github.com/Deltares/Wflow.jl/blob/master/test/sbm_gwf_config.toml)):

```toml
[model]
drains = true
```

### Recharge
The recharge flux ``Q_{r}`` to the aquifer is calculated as follows:

```math
Q_{r} = R \, A
```
with ``R`` the recharge rate [L T``^{-1}``] and ``A`` the area [L``^2`` ] of the aquifer
cell.

The Table in the Groundwater flow [recharge boundary condition](@ref gwf_recharge_params)
section of the Model parameters section provides the parameters of the struct `Recharge`.
Parameters that can be set directly from the static input data (netCDF) are marked in this
Table.

The recharge flux ``Q_r`` is an output variable (field `flux` of struct `Recharge`), and is
used to update the total flux in a cell where recharge occurs. For the model `SBM +
Groundwater flow`, the recharge rate from the vertical SBM concept `recharge` [mm] is used
to update the `rate` field of the `Recharge` struct each time step. The `rate` field is
multiplied by the `area` field of the aquifer.   

### Head boundary
This boundary is a fixed head with time (not affected by the model stresses over time))
outside of the model domain, and is generally used to avoid an unnecessary extension of the
model domain to the location of the fixed boundary (for example a large lake). The flux from
the boundary ``Q_{hb}`` [L``^3`` T``^{-1}``] is calculated as follows:

```math
Q_{hb} = C_{hb} (\phi_{hb} - \phi)
```
with ``C_{hb}`` the conductance of the head boundary [L``^2`` T``^{-1}``], ``\phi_{hb}`` the
head [L] of the head boundary and  ``\phi`` the head of the aquifer cell.

The Table in the Groundwater flow [head boundary condition](@ref gwf_headboundary_params)
section of the Model parameters provides the parameters of the struct `HeadBoundary`.

The head boundary flux ``Q_{hb}`` is an output variable (field `flux` of struct
`HeadBoundary`), and is used to update the total flux in a cell where this type of boundary
occurs. The parameter Head ``\phi_{hb}`` can be specified as a fixed or time dependent
value.

!!! note 
    This boundary is not (yet) part of the model `SBM + Groundwater flow`.

### Well boundary
A volumetric well rate [L``^3`` T``^{-1}``] can be specified as a boundary condition.

The Table in the [well boundary condition](@ref well_boundary_params) section of the Model
parameters provides the parameters of the struct `Well`.

The volumetric well rate ``Q_{well}`` can be can be specified as a fixed or time dependent
value. If a cell is dry, the actual well flux `flux` is set to zero (see also the last note
on this page).

!!! note 
    This boundary is not (yet) part of the model `SBM + Groundwater flow`.

!!! note 
    For an unconfined aquifer the boundary fluxes are checked, in case of a dry aquifer cell
    a negative flux is not allowed.

## References
+ Chu, W. S., & Willis, R. (1984). An explicit finite difference model for unconfined
  aquifers. Groundwater, 22(6), 728-734.
+ Goode, D. J., & Appel, C. A. (1992). Finite-Difference Interblock Transmissivity for
  Unconﬁned Aquifers and for Aquifers having Smoothly Varying Transmissivity Water-resources
  investigations report, 92, 4124. 
+ Johnson, A. I. (1967), Specific yield: compilation of specific yields for various
  materials, Water Supply Paper 1662-D, Washington, D.C.: U.S. Government Printing Office,
  p. 74, doi:10.3133/wsp1662D.
+ Langevin, C.D., Hughes, J.D., Banta, E.R., Niswonger, R.G., Panday, Sorab, and Provost,
  A.M., 2017, Documentation for the MODFLOW 6 Groundwater Flow Model: U.S. Geological Survey
  Techniques and Methods, book 6, chap. A55, 197 p., https://doi.org/10.3133/tm6A55.
+ McDonald, M.G., Harbaugh, A.W., Orr, B.R., and Ackerman, D.J., 1991, A method of
  converting no-flow cells to variable-head cells for the U.S. Geological Survey modular
  finite-difference groundwater flow model: U.S. Geological Survey Open-File Report 91-536,
  99 p.
+ Niswonger, R.G., Panday, Sorab, and Ibaraki, Motomu, 2011, MODFLOW-NWT, A Newton
  formulation for MODFLOW-2005: U.S. Geological Survey Techniques and Methods 6-A37, 44 p.