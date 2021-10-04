## [River flow](@id local-inertial_river_params)
The Table below shows the parameters (fields) of struct `ShallowWaterRiver`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[input.lateral.river]`, to map the
internal model parameter to the external netCDF variable.

|  parameter  | description  	  | unit  | default |
|:--------------- | ------------------| ----- | -------- |
| **`mannings_n`**    |  Manning's roughness at edge/link| s m``^{-\frac{1}{3}}`` | 0.036 |
| **`width`**    |  river width | m | - |
| **`z`**    |  river bed elevation | m | - |
| **`length`**    |  river length | m | - |
| `n`   |  number of cells | - | - |
| `ne`    |  number of edges/links | - | - |
| `g`    |  acceleration due to gravity | m s``^{-2}`` | - |
| `α`    |  stability coefficient (Bates et al., 2010) | - | - |
| `h_thresh`    |  depth threshold for calculating flow | m | - |
| `Δt`    |  model time step | s | - |
| `q`    |  river discharge (subgrid channel) | m``^3`` s``^{-1}`` | - |
| `q_av`    |  average river discharge (subgrid channel) | m``^3`` s``^{-1}`` | - |
| `zmax`    | maximum channel bed elevation | m | - |
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
| `inwater0`    | lateral inflow at previous time step | m``^3`` s``^{-1}`` | - |
| `bankvolume`    | bank volume | m``^3`` | - |
| `bankheight`    | bank height | m | - |
| `froude_limit`    | if true a check is performed if froude number > 1.0 (algorithm is modified) | - | - |
| `reservoir_index`   |  map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field) | - | - |
| `lake_index`   |  map cell to 0 (no lake) or i (pick lake i in lake field) | - | - |
| `reservoir`    | an array of reservoir models `SimpleReservoir` | - | - |
| `lake` | an array of lake models `NaturalLake` | - | - |