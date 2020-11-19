# Kinematic wave

## Surface routing
The main flow routing scheme available in Wflow.jl is the kinematic wave approach for
channel and overland flow, assuming that the topography controls water flow mostly. The
kinemative wave equations are (Chow, 1988): ``\dfrac{dQ}{dx} + \dfrac{dA}{dt} = q`` and ``A
= \alpha Q^{\beta}``. These equations can then be combined as a function of streamflow only:
```math
    \dfrac{dQ}{dx} + \alpha \beta Q^{\beta - 1} \dfrac{dQ}{dt} = q
```
where ``Q`` is the surface runoff in the kinematic wave [m``^3``/s], ``x`` is the length of
the runoff pathway [m], ``A`` is the cross-section area of the runoff pathway [m``^{2}``],
``t`` is the integration timestep [s] and ``alpha`` and ``\beta`` are coefficients.

These equations are solved with a nonlinear scheme using Newton’s method and can also be
iterated depending on the  model space and time resolution. By default, the iterations are
performed until a stable solution is reached (``\epsilon < 10^{-12}``). For larger models,
the number of iterations can also be fixed for to a specific sub-timestep (in seconds) for
both overland and channel flows to improve simulation time. To enable (fixed or not)
iterations of the kinematic wave the following lines can be inserted in the TOML file of the
model:

```toml
[model]
# Enable iterations of the kinematic wave
kin_wave_iteration = true 
```
## Subsurface flow routing
In the SBM model the kinematic wave approach is used to route subsurface flow laterally. The
saturated store ``S`` can be drained laterally by saturated downslope subsurface flow per
unit width of slope ``w`` [mm] according to:
```math
    q=\frac{K_{0}\mathit{tan(\beta)}}{f}(e^{(-fz_{i})}-e^{(-fz_{t})})
```
where ``\beta`` is element slope angle [deg.], ``q`` is subsurface flow [mm``^{2}``/t],
``K_{0}`` is the saturated hydraulic conductivity at the soil surface [mm/t], ``z_{i}`` is
the water table depth [mm],``z_{t}`` is total soil depth [mm], and ``f`` is a scaling
parameter [mm``^{-1}``]:
```math
    f=\frac{\theta_{s}-\theta_{r}}{M},\,
```
where ``\theta_{s}`` is saturated water content [mm/mm] and ``\theta_{r}`` is residual water
content [mm/mm] and ``M`` represents a model parameter [mm], that determines the decrease of
vertical saturated conductivity with depth.

Combining with the following continuity equation:
```math
    (\theta_s-\theta_r)\frac{\partial h}{\partial t} = -w\frac{\partial q}{\partial x} + wr
```
where ``h`` is the water table height [mm], ``x`` is the distance downslope [mm], and ``r``
is the netto input rate [mm/t] to the saturated store. Substituting for ``h (\frac{\partial
q}{\partial h})``, gives:
```math 
  w \frac{\partial q}{\partial t} = -cw\frac{\partial q}{\partial x} + cwr
```

where celerity ``c = \frac{K_{0}\mathit{tan(\beta)}}{(\theta_s-\theta_r)} e^{(-fz_{i})}``

The kinematic wave equation for lateral subsurface flow is solved iteratively using Newton's
method.

!!! note For the lateral subsurface flow kinematic wave the model timestep is not adjusted.
    For certain model timestep and model grid size combinations this may result in loss of
    accuracy.

## Subcatchment flow
Normally the the kinematic wave is continuous throughout the model. By using the the `pits`
entry in the model and input sections of the TOML file all flow is at the subcatchment only
(upstream of the pit locations, defined by the netCDF variable `wflow_pits` in the example
below) and no flow is transferred from one subcatchment to another. This can be convenient
when connecting the result of the model to a water allocation model such as Ribasim.

```toml
[input]
# these are not directly part of the model
pits = "wflow_pits"

[model]
pits = true
```

## Limitations
The kinematic wave approach for channel, overland and lateral subsurface flow, assumes that
the topography controls water flow mostly. This assumption holds for steep terrain, but in
less steep terrain the hydraulic gradient is likely not equal to the surface slope
(subsurface flow), or pressure differences and inertial momentum cannot be neglected
(channel and overland flow). In addition, while the kinemative wave equations are solved
with a nonlinear scheme using Newton's method (Chow, 1988), other model equations are solved
through a simple explicit scheme. In summary the following limitations apply:

+ Channel flow, and to a lesser degree overland flow, may be unrealistic in terrain that is
  not steep, and where pressure forces and inertial momentum are important.

+ The lateral movement of subsurface flow may be very wrong in terrain that is not steep.
