# [Kinematic wave] (@id kin_wave)

## Surface routing
The main flow routing scheme available in Wflow.jl is the kinematic wave approach for
channel and overland flow, assuming that the topography controls water flow mostly. The
kinematic wave equations are (Chow, 1988):
```math
  \dfrac{dQ}{dx} + \dfrac{dA}{dt} = q \\~\\
   A = \alpha Q^{\beta}
```
These equations can then be combined as a function of streamflow only:
```math
    \dfrac{dQ}{dx} + \alpha \beta Q^{\beta - 1} \dfrac{dQ}{dt} = q
```
where ``Q`` is the surface runoff in the kinematic wave [m``^3``/s], ``x`` is the length of
the runoff pathway [m], ``A`` is the cross-section area of the runoff pathway [m``^{2}``],
``t`` is the integration timestep [s] and ``\alpha`` and ``\beta`` are coefficients.

These equations are solved with a nonlinear scheme using Newton's method and can also be
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
# Fixed sub-timestep for iterations of channel flow (river cells)
kw_river_tstep = 900
# Fixed sub-timestep for iterations of overland flow (land cells)
kw_land_tstep = 3600
```

The ``\alpha`` parameter of the kinematic wave is fixed. To estimate the wetted perimeter
for the calculation of the ``\alpha`` parameter a bankfull river depth map (default value
is 1.0 m) for the river can be provided as follows:

```toml
[input.lateral.river]
bankfull_depth = "wflow_riverdepth"
```

The wetted perimeter of the river is based on half bankfull river depth. For the land part the
wetted perimeter is based on the flow width.

Simplified [reservoir and lake](@ref reservoir_lake) models can be included as part of the
river kinematic wave network.

## Inflow
External water (supply/abstraction) `inflow` [m``^3`` s``^{-1}``]  can be added to the
kinematic wave for surface water routing, as a cyclic parameter or as part of forcing (see
also [Input section](@ref)).

## Subsurface flow routing
In the SBM model the kinematic wave approach is used to route subsurface flow laterally.
Different vertical hydraulic conductivity depth profiles are possible as part of the
vertical [SBM](@ref soil) concept, and these profiles (after unit conversion) are also used
to compute lateral subsurface flow. The following profiles (see [SBM](@ref soil) for a
detailed description) are available:
- `exponential` (default)
- `exponential_constant`
- `layered`
- `layered_exponential` 
For the profiles `exponential` and `exponential_constant`, the saturated store ``S`` is
drained laterally by saturated downslope subsurface flow for a slope with width ``w`` [m]
according to:
```math
    Q = \begin{cases}
    \frac{K_0\tan(\beta)}{f}\left(e^{(-fz_{i})}-e^{(-fz_\mathrm{exp})}\right) w + 
    K_0e^{(-fz_\mathrm{exp})}(z_t-z_\mathrm{exp})\tan(\beta) w & \text{if $z_i < z_\mathrm{exp}$}\\
    \\
    K_0e^{(-fz_\mathrm{exp})}(z_t - z_i)\tan(\beta) w & \text{if $z_i \ge z_\mathrm{exp}$},
    \end{cases}
```
where ``\beta`` is element slope angle, ``Q`` is subsurface flow [m``^{3}`` d``^{-1}``],
``K_0`` is the saturated hydraulic conductivity at the soil surface [m d``^{-1}``], ``z_i``
is the water table depth [m], ``z_{t}`` is the total soil depth [m], ``f`` is a scaling
parameter [m``^{-1}``] that controls the decrease of ``K_0`` with depth and
``z_\mathrm{exp}`` [m] is the depth from soil surface for which the exponential decline of
``K_0`` is valid. For the `exponential` profile, ``z_\mathrm{exp}`` is equal to ``z_t``.

Combining with the following continuity equation:
```math
    (\theta_s-\theta_r)w\frac{\partial h}{\partial t} = -\frac{\partial Q}{\partial x} + wr
```
where ``h`` is the water table height [m], ``x`` is the distance downslope [m], and ``r`` is
the net input rate [m d``^{-1}``] to the saturated store. Substituting for ``h
(\frac{\partial Q}{\partial h})``, gives:
```math
  \frac{\partial Q}{\partial t} = -c\frac{\partial Q}{\partial x} + cwr
```

where celerity ``c`` is calculated as follows:
```math
    c = \begin{cases}
    \frac{K_0e^{(-fz_{i})}\tan(\beta)}{(\theta_s-\theta_r)} 
    + \frac{K_0e^{(-fz_\mathrm{exp})}\tan(\beta)}{(\theta_s-\theta_r)}  & \text{if $z_i < z_\mathrm{exp}$}\\
    \\
    \frac{K_0e^{(-fz_\mathrm{exp})}\tan(\beta)}{(\theta_s-\theta_r)} & \text{if $z_i \ge z_\mathrm{exp}$}.
    \end{cases}
```

For the `layered` and `layered_exponential` profiles the equivalent horizontal hydraulic
conductivity ``K_h`` [m d``^{-1}``] is calculated for water table height ``h = z_t-z_i``
[m], and lateral subsurface flow is calculated as follows:
```math
  Q = K_h h \tan(\beta) w,
```
and celerity ``c`` is given by:
```math
    c = \frac{K_h \tan(\beta)}{(\theta_s-\theta_r)}.
```

The kinematic wave equation for lateral subsurface flow is solved iteratively using Newton's
method.

!!! note
    For the lateral subsurface flow kinematic wave the model timestep is not adjusted.
    For certain model timestep and model grid size combinations this may result in loss of
    accuracy.

## Multi-Threading
The kinematic wave calculations for surface - and subsurface flow routing can be executed in
parallel using multiple threads. In the model section of the TOML file, a minimum stream
order can be provided to define subbasins for the river (default is 6) and land domain
(default is 5). Subbasins are created at all confluences where each branch has a minimal
stream order. Based on the subbasins a directed acyclic graph is created that controls the
order of execution and which subbasins can run in parallel.

```toml
[model]
min_streamorder_river = 5 # minimum stream order to delineate subbasins for river domain, default is 6
min_streamorder_land = 4  # minimum stream order to delineate subbasins for land domain, default is 5
```

## Subcatchment flow
Normally the the kinematic wave is continuous throughout the model. By using the `pits`
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
(channel and overland flow). In addition, while the kinematic wave equations are solved
with a nonlinear scheme using Newton's method (Chow, 1988), other model equations are solved
through a simple explicit scheme. In summary the following limitations apply:

+ Channel flow, and to a lesser degree overland flow, may be unrealistic in terrain that is
  not steep, and where pressure forces and inertial momentum are important.

+ The lateral movement of subsurface flow may be very wrong in terrain that is not steep.

## External inflows
External inflows, for example water supply or abstractions, can be added to the kinematic
wave via the `inflow` variable. For this, the user can supply a 2D map of the inflow, as a
cyclic parameter or as part of forcing (see also [Input section](@ref)). These inflows are
added or abstracted from the upstream inflow `qin` before running the kinematic wave to
solve the impact on resulting `q`. In case of a negative inflow (abstractions), a minimum of
zero is applied to the upstream flow `qin`.

## References
+ Chow, V., Maidment, D. and Mays, L., 1988, Applied Hydrology. McGraw-Hill Book Company,
  New York.
