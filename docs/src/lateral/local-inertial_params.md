## [River flow](@id local-inertial_river_params)
The Table below shows the parameters (fields) of struct `ShallowWaterRiver`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[input.lateral.river]`, to map the
internal model parameter to the external netCDF variable. The parameter river bed elevation
`zb` is based on the bankfull elevation and depth input data:

```toml
[input.lateral.river]
bankfull_elevation = "RiverZ"
bankfull_depth = "RiverDepth"
```

|  parameter  | description  	  | unit  | default |
|:--------------- | ------------------| ----- | -------- |
| **`mannings_n`**    |  Manning's roughness at edge/link| s m``^{-\frac{1}{3}}`` | 0.036 |
| **`width`**    |  river width | m | - |
| `zb`    |  river bed elevation | m | - |
| **`length`**    |  river length | m | - |
| `n`   |  number of cells | - | - |
| `ne`    |  number of edges/links | - | - |
| `g`    |  acceleration due to gravity | m s``^{-2}`` | - |
| `α`    |  stability coefficient (Bates et al., 2010) | - | 0.7 |
| `h_thresh`    |  depth threshold for calculating flow | m | 0.001 |
| `Δt`    |  model time step | s | - |
| `q`    |  river discharge (subgrid channel) | m``^3`` s``^{-1}`` | - |
| `q_av`    |  average river discharge (subgrid channel) | m``^3`` s``^{-1}`` | - |
| `zb_max`    | maximum channel bed elevation | m | - |
| `h`    | water depth | m | - |
| `η_max`    | maximum water elevation | m | - |
| `hf`    | water depth at edge/link | m | - |
| `h_av`    | average water depth | m | - |
| `length_at_link`    | river length at edge/link | m | - |
| `width_at_link`    | river width at edge/link | m | - |
| `a`    | flow area at edge/link | m``^2`` | - |
| `r`    | wetted perimeter at edge/link | m | - |
| `volume`    | river volume | m``^3`` | - |
| `error`    | error volume | m``^3`` | - |
| `inwater`    | lateral inflow | m``^3`` s``^{-1}`` | - |
| `inflow`        | external inflow (abstraction/supply/demand) | m``^3`` s``^{-1}``| 0.0 |
| `bankfull_volume`    | bankfull volume | m``^3`` | - |
| **`bankfull_depth`**    | bankfull depth | m | - |
| `froude_limit`    | if true a check is performed if froude number > 1.0 (algorithm is modified) | - | - |
| `reservoir_index`   |  map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field) | - | - |
| `lake_index`   |  map cell to 0 (no lake) or i (pick lake i in lake field) | - | - |
| `reservoir`    | an array of reservoir models `SimpleReservoir` | - | - |
| `lake` | an array of lake models `NaturalLake` | - | - |

## [Overland flow](@id local-inertial_land_params)
The Table below shows the parameters (fields) of struct `ShallowWaterLand`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[input.lateral.land]`, to map the
internal model parameter to the external netCDF variable.

|  parameter  | description  	  | unit  | default |
|:--------------- | ------------------| ----- | -------- |
| `n`    |  number of cells | - | - |
| `xl`| cell length x direction | m | - | 
| `yl`| cell length y direction | m | - | 
| `xwidth`| effective flow width x direction (floodplain) | m | - | 
| `ywidth`| effective flow width y direction (floodplain) | m | - |
| `g` | acceleration due to gravity | m s``^{-2}`` | - |
| `θ` | weighting factor (de Almeida et al., 2012) | - | 0.8 |
| `α`    |  stability coefficient (Bates et al., 2010) | - | 0.7 |
| `h_thresh`    |  depth threshold for calculating flow | m | 0.001 |
| `Δt`   |  model time step| s | - |
| `qy0`  |  flow in y direction at previous time step| m``^3`` s``^{-1}`` | - |
| `qx0`  |  flow in x direction at previous time step| m``^3`` s``^{-1}`` | - |
| `qx`  |  flow in x direction | m``^3`` s``^{-1}`` | - |
| `qy`  |  flow in y direction | m``^3`` s``^{-1}`` | - |
| `zx_max`  |  maximum cell elevation (x direction) | m | - |
| `zy_max`  |  maximum cell elevation (y direction) | m | - |
| **`mannings_n`**  |  Manning's roughness | s m``^{-\frac{1}{3}}`` | 0.072 |
| `volume`  |  total volume of cell | m``^3`` | - |
| `error`  |  error volume | m``^3`` | - |
| `runoff`  |  runoff from hydrological model | m``^3`` s``^{-1}`` | - |
| `h`  |  water depth of cell | m``^3`` s``^{-1}`` | - |
| **`z`**  |  elevation of cell | m | - |
| `froude_limit`  |  if true a check is performed if froude number > 1.0 (algorithm is modified)| - | - |
| `rivercells`  |  river cells| - | - |
| `h_av`  | average water depth| m | - |