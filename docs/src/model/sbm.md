# wflow\_sbm

Wflow\_sbm represents a family of hydrological models derived from the CQflow model (Köhler
et al.,2006) that has been applied in various countries, most notably in Central America.
The models have the vertical SBM concept in common. The soil part of wflow\_sbm is largely
based on the Topog\_SBM model but has had considerable changes over time. Topog\_SBM is
specifically designed to simulate fast runoff processes in small catchments while wflow\_sbm
model can be applied more widely. The main differences are for the vertical concept SBM of
wflow\_sbm:

- The unsaturated zone can be split-up in different layers
- The addition of evapotranspiration losses
- The addition of a capillary rise

The vertical SBM concept is explained in more detail in the following section [SBM vertical
concept](@ref).

Topog\_SBM uses an element network based on contour lines and trajectories for water
routing. Wflow\_sbm models differ in how the lateral components:
- river
- land
- subsurface  

are solved.

## SBM + Kinematic wave
For the lateral components of this wflow\_sbm model water is routed over a D8 network, and
the kinematic wave approach is used for river, overland and lateral subsurface flow. This is
described in more detail in the section [Kinematic wave](@ref).

Overview of the different processes and fluxes in the wflow_sbm model:

![wflow_sbm model](../images/wflow_sbm_soil.png)

Below the mapping for this wflow\_sbm model (type `sbm`) to the vertical SBM concept
(instance of `struct SBM`) and the different lateral concepts is presented. For an
explanation about the type parameters between curly braces after the `struct` name see
[Vertical and lateral models](@ref).

```julia
vertical => struct SBM{T,N,M}
lateral.subsurface => struct LateralSSF{T}
lateral.land => struct SurfaceFlow{T,R,L}
lateral.river => struct SurfaceFlow{T,R,L}
lateral.river.lake => struct NaturalLake{T} # optional
lateral.river.reservoir => struct SimpleReservoir{T} # optional
```

## SBM + Local inertial river
By default the [SBM + Kinematic wave](@ref) model (type `sbm`) uses the kinematic wave
approach for river flow. There is also the option to use the local inertial model for
river flow, by providing the following in the TOML file:  

```toml
[model]
river_routing = "local-inertial"
```

Only the mapping for the river component changes, as shown below. For an explanation about
the type parameters between curly braces after the `struct` name see [Vertical and lateral
models](@ref). 

```julia
lateral.river => struct ShallowWaterRiver{T,R,L}
```
The local inertial approach is described in more detail in the section [Local inertial
model](@ref).

## SBM + Groundwater flow
For river and overland flow the kinematic wave approach over a D8 network is used for this
wflow\_sbm model. For the subsurface domain, an unconfined aquifer with groundwater flow in
four directions (adjacent cells) is used. This is described in more detail in the section
[Groundwater flow](@ref lateral_gwf).

Below the mapping for this wflow\_sbm model (type `sbm_gwf`) to the vertical SBM concept
(instance of `struct SBM`) and the different lateral concepts. For an explanation about the
type parameters between curly braces after the `struct` name see [Vertical and lateral
models](@ref).

```julia
vertical => struct SBM{T,N,M}
lateral.subsurface.flow => struct GroundwaterFlow{A, B}
lateral.subsurface.recharge => struct Recharge{T} <: AquiferBoundaryCondition
lateral.subsurface.river => struct River{T} <: AquiferBoundaryCondition
lateral.subsurface.drain => struct Drainage{T} <: AquiferBoundaryCondition # optional
lateral.land => struct SurfaceFlow{T,R,L}
lateral.river => struct SurfaceFlow{T,R,L}
lateral.river.lake => struct NaturalLake{T} # optional
lateral.river.reservoir => struct SimpleReservoir{T} # optional
```

## References
+ Köhler, L., Mulligan, M., Schellekens, J., Schmid, S., Tobón, C., 2006, Hydrological
  impacts of converting tropical montane cloud forest to pasture, with initial reference to
  northern Costa Rica. Final Technical Report DFID‐FRP Project No. R799.
