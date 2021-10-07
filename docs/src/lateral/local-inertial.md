# Local inertial model

## River routing
The local inertial approximation of shallow water flow neglects only the convective
acceleration term in the Saint-Venant momentum conservation equation. The numerical solution
of the local inertial approximation on a staggered grid is as follows (Bates et al., 2010):

```math
Q_{t+\Delta t} = \frac{Q_t - g A_t \Delta t S_t}{(1+g\Delta t n^2 |Q_t| / (R_t^{4/3} A_t))}
``` 
where ``Q_{t+\Delta t}`` is the river flow [m``^3``/s] at time step ``t+\Delta t``, ``g`` is
acceleration due to gravity [m/s``^2``], ``A_t`` is the cross sectional flow area at the
previous time step, ``R_t`` is the hydraulic radius at the previous time step, ``Q_t`` is
the river flow [m``^3``/s] at the previous time step, ``S_t`` is the water surface slope at
the previous time step and ``n`` is the Manning's roughness coefficient [m``^{-1/3}`` s].

The momentum equation is applied to each link between two river grid cells, while the
continuity equation over ``\Delta t`` is applied to each river cell:

```math
h^{t+\Delta t} = h^t + \Delta t \frac{Q^{t+\Delta t}_{src} - Q^{t+\Delta t}_{dst}}{A}
``` 
where ``h^{t+\Delta t}`` is the water depth [m] at time step ``t+\Delta t``, ``h^t`` is the
water depth [m] at the previous time step, ``A`` is the river area [m``^2``] and ``Q_{src}``
and ``Q_{dst}`` represent river flow [m``^3``/s] at the upstream and downstream link of the
river cell, respectively.

The model time step ``\Delta t`` for the local inertal model is estimated based on the
Courant-Friedrichs-Lewy condition (Bates et al., 2010):

```math
\Delta t = min(\alpha \frac{\Delta x_i}{\sqrt{(gh_i)}})
```

where ``\sqrt{(gh_i)}`` is the wave celerity for river cell ``i`` , ``\Delta x_i`` is the
river length [m] for river cell ``i`` and ``\alpha`` is a coefficient (typically between 0.2
and 0.7) to enhance the stability of the simulation.

In the TOML file the following properties related to the local inertial model can be
provided for the `sbm` model type:

```toml
[model]
river_routing = "local-inertial" # default is kinematic-wave
inertial_flow_alpha = 0.5      # alpha coefficient for model stability (default = 0.7)
froude_limit = true            # default is true, limit flow to subcritical-critical according to Froude number
h_thresh = 0.1                 # water depth [m] threshold for calculating flow between cells (default = 1e-03)
riverlength_bc = 1000          # river length [m] for boundary points (default = 1e05)
```

The momentum equation is most stable for low slope enviroments, and to keep the simulation
stable for (partly) steep environments the `froude_limit` option can be set to true. This
setting limits flow to subcritical-critical conditions based on the Froude number ($\le 1$),
similar to Coulthard et al. (2013) in the CAESAR-LISFLOOD model and Adams et al. (2017) in
the Landlab v1.0 OverlandFlow component.

For the downstream boundary condition (ghost point) a fixed water surface elevation of 0 m
is used. The river length [m] of the boundary cell can be set through the TOML file with
`riverlength_bc`, the default value is 100 km. 

Simplified reservoirs and lakes models can be included as part of the local inertial model, 
see also [Reservoirs and Lakes](@ref).

## Multi-Threading
The local inertial model for river flow can be executed in parallel using multiple threads.

## References
+ Adams, J. M., Gasparini, N. M., Hobley, D. E. J., Tucker, G. E., Hutton, E. W. H.,
  Nudurupati, S. S., and Istanbulluoglu, E., 2017, The Landlab v1.0 OverlandFlow component:
  a Python tool for computing shallow-water flow across watersheds, Geosci. Model Dev., 10,
  1645–1663, <https://doi.org/10.5194/gmd-10-1645-2017>. 
+ Bates, P. D., M. S. Horritt, and T. J. Fewtrell, 2010, A simple inertial formulation of
  the shallow water equations for efficient two-dimensional flood inundation modelling, J.
  Hydrol., 387, 33–45, doi:10.1016/ j.jhydrol.2010.03.027.
+ Coulthard, T. J., Neal, J. C., Bates, P. D., Ramirez, J., de Almeida, G. A. M., and
  Hancock, G. R., 2013, Integrating the LISFLOOD-FP 2- D hydrodynamic model with the CAESAR
  model: implications for modelling landscape evolution, Earth Surf. Proc. Land., 38,
  1897–1906, doi:10.1002/esp.3478.
