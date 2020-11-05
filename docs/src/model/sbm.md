# wflow\_sbm

Wflow\_sbm represents a family of hydrological models that derived from the
CQflow model (Köhler et al.,2006) that has been applied in various countries,
most notably in Central America. The models have the vertical SBM concept in
common. The soil part of wflow\_sbm is largely based on the Topog\_SBM model but
has had considerable changes over time. Topog\_SBM is specifically designed to
simulate fast runoff processes in small catchments while wflow\_sbm model can be
applied more widely. The main differences are for the vertical concept SBM of
wflow\_sbm:

- The unsaturated zone can be split-up in different layers
- The addition of evapotranspiration losses
- The addition of a capillary rise

The vertical SBM concept is explained in more detail in [SBM vertical
concept](@ref).

Topog\_SBM uses an element network based on contour lines and trajectories for
water routing. Wflow\_sbm models differ in how the lateral components:
- river
- land
- subsurface  

are solved.

## SBM + Kinematic wave
For the lateral components of this wflow\_sbm model water is routed over a D8
network, and the kinematic wave approach is used for river, overland and lateral
subsurface flow. This is described in more detail [Kinematic wave](@ref).

Overview of the different processes and fluxes in the wflow_sbm
model:

![wflow_sbm model](../images/wflow_sbm_soil.png)

## SBM + Groundwater flow
For river and overland flow the kinematic wave approach over a D8 network is
used for this wflow\_sbm model. For the subsurface domain, an unconfined aquifer
with groundwater flow in four directions (adjacent cells) is used. This is
described in more detail [Groundwater flow](@ref).
