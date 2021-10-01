## Surface flow
The Table below shows the parameters (fields) of struct `SurfaceFlow`, including a
description of these parameters, the unit, and default value if applicable. `SurfaceFlow` is
used for river and overland flow. The parameters in bold represent model parameters that can
be set through static input data (netCDF), and can be listed in the TOML configuration file
under `[input.lateral.river]` and `[input.lateral.land]`, for river and overland flow
respectively, to map the internal model parameter to the external netCDF variable. For river
flow three additionally parameters can be set, `dl` (river length), `width` (river width)
and bankfull water height `h_bankfull` as follows, through the TOML file:

```toml
[input.lateral.river]
length = "wflow_riverlength"
width = "wflow_riverwidth"
h_bankfull = "river_bankfullheight"
```
[^1]: default value for Manning's roughness `n`, river = 0.036; land = 0.072
[^2]: only applicable for river domain

|  parameter  | description  	  | unit  | default |
|:--------------- | ------------------| ----- | -------- |
| `β`             |  constant in Manning's equation | - | - |
| `dl`            |  length      | m     |  -      |
| **`n`**             | Manning's roughness | s m``^{\frac{1}{3}}``| [^1] |
| **`sl`**         |  slope       | m m``^{-1}``| - |
| `width`         |  width       | m          | - |
| `q`             | discharge     | m``^3`` s``^{-1}``| - |
| `qin`           | inflow from upstream cells | m``^3`` s``^{-1}``| - |
| `q_av`          | average discharge | m``^3`` s``^{-1}``| - |
| `qlat`          | lateral inflow per unit length  | m``^2`` s``^{-1}``| - |
| `inwater`       | lateral inflow | m``^3`` s``^{-1}``| - |
| `inflow`        | external inflow (abstraction/supply/demand) | m``^3`` s``^{-1}``| 0.0 |
| `volume`        | kinematic wave volume |m``^3``| - |
| `h`             | water level | m | - |
| `h_av`             | average water level | m | - |
| `h_bankfull`     | bankfull water level  | m | 1.0 | 
| `Δt`             | model time step | s | - |
| `its`             | number of fixed iterations | - | - |
| `alpha_pow`             | used in the power part of ``\alpha`` | - | - |
| `alpha_term`             | term used in computation of ``\alpha`` | - | - |
| `α`            | constant in momentum equation ``A = \alpha Q^{\beta}`` | s``^{\frac{3}{5}}`` m``^{\frac{1}{5}}`` | - |
| `cel`            | celerity of kinematic wave | m s``^{-1}`` | - |
| `to_river`            | part of overland flow that flows to the river | m s``^3`` | - |
| `rivercells`            | location of river cells (0 or 1) | - | - |
| `wb_pit`            |  location (0 or 1) of a waterbody (wb, reservoir or lake) | - | - |
| `reservoir_index`   |  map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field) | - | - |
| `lake_index`   |  map cell to 0 (no lake) or i (pick lake i in lake field) | - | - |
| `reservoir` [^2]    | an array of reservoir models `SimpleReservoir` | - | - |
| `lake` [^2]    | an array of lake models `NaturalLake` | - | - |
| `kinwave_it`   | boolean for kinematic wave iterations | - | false |

## [Reservoirs](@id reservoir_params)
The Table below shows the parameters (fields) of struct `SimpleReservoir`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[input.lateral.river.reservoir]`, to map
the internal model parameter to the external netCDF variable.

Two parameters reservoir coverage `areas` and the outlet of reservoirs (unique id) `locs`
that are not part of the `SimpleReservoir` struct are also required, and can be set as
follows through the TOML file: 

```toml
[input.lateral.river.reservoir]
areas = "wflow_reservoirareas"
locs = "wflow_reservoirlocs"
```

|  parameter | description  	        | unit | default |
|:---------------| --------------- | ---------------------- | ----- |    
| **`area`**        | area | m``^2`` | - |
| **`demandrelease`** ( `demand`) | minimum (environmental) flow released from reservoir  | m``^3`` s``^{-1}``| - |
| **`maxrelease`** | maximum amount that can be released if below spillway | m``^3`` s``^{-1}`` | - |
| **`maxvolume`** | maximum storage (above which water is spilled) | m``^3`` | - |
| **`targetfullfrac`** | target fraction full (of max storage)| -        | - |
| **`targetminfrac`** | target minimum full fraction (of max storage) | -  | - |
| `Δt`             | model time step     | s | - |
| `volume`             | volume     | m``^3`` | - |
| `inflow`             | total inflow into reservoir | m``^3`` | - |
| `outflow`             | outflow into reservoir | m``^3`` s``^{-1}`` | - |
| `totaloutflow`        | total outflow into reservoir | m``^3`` | - |
| `percfull`             | fraction full (of max storage) | - | - |
| `precipitation`             | outflow into reservoir | mm Δt⁻¹ | - |
| `evaporation`             | outflow into reservoir | mm Δt⁻¹ | - |

## [Lakes](@id lake_params)
The Table below shows the parameters (fields) of struct `NaturalLake`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[input.lateral.river.lake]`, to map the
internal model parameter to the external netCDF variable.

Two parameters lake coverage `areas` and the outlet of lakes (unique id) `locs` that are not
part of the `NaturalLake` struct are also required, and can be set as follows through the
TOML file:

```toml
[input.lateral.river.lake]
areas = "wflow_lakeareas"
locs = "wflow_lakelocs"
```

|  parameter | description  	        | unit | default |
|:---------------| --------------- | ---------------------- | ----- |    
| **`area`**         | area| m``^2`` | - |
| **`b`**       | Rating curve coefficient  | - | - |
| **`e`**      | Rating curve exponent | -  | - |
| **`outflowfunc`**     | type of lake rating curve | - | - |
| **`storfunc`** | type of lake storage curve| - | - |
| **`threshold`**  | water level threshold ``H_0`` below that level outflow is zero | m  | - |
| **`waterlevel`**  | waterlevel ``H`` of lake | m  | - |
| **`lowerlake_ind`** (`linkedlakelocs`) | Index of lower lake (linked lakes) | - | 0 |
| **`sh`** | data for storage curve | - | - |
| **`hq`** | data rating curve | - | - |
| `Δt`             | model time step     | s |  - |
| `inflow` | total inflow to the lake | m``^3``  | - |
| `storage` | storage lake | m``^3``  | - |
| `outflow` | outflow lake | m``^3`` s``^{-1}``  | - |
| `totaloutflow` | total outflow lake | m``^3``  | - |
| `precipitation` | average precipitation for lake area | mm Δt⁻¹  | - |
| `evaporation` | average precipitation for lake area | mm Δt⁻¹  | - |

## Lateral subsurface flow
The Table below shows the parameters (fields) of struct `LateralSSF`, including a
description of these parameters, the unit, and default value if applicable. 

| parameter | description  	        | unit | default |
|:---------------| --------------- | ---------------------- | ----- |    
| `kh₀`          | horizontal hydraulic conductivity at soil surface  | m Δt``^{-1}`` |  - |
| `f` | a scaling parameter (controls exponential decline of kh₀) | m``^{-1}``  | - |
| `soilthickness` | soil thickness | m  | - |
| `θₛ` | saturated water content (porosity) | -  | - |
| `θᵣ` | residual water content  | -  | - |
| `t` | time step | Δt s  | - |
| `Δt` | model time step | s  | - |
| `βₗ` | slope | -  | - |
| `dl` | drain length | m | - |
| `dw` | drain width | m | - |
| `zi` | pseudo-water table depth (top of the saturated zone) | m | - |
| `exfiltwater` | exfiltration (groundwater above surface level, saturated excess conditions) | m | - |
| `recharge` | net recharge to saturated store  | m | - |
| `ssf` | subsurface flow  | m``^3`` Δt``{-1}``  | - |
| `ssfin` | inflow from upstream cells | m``^3`` Δt``{-1}``  | - |
| `ssfmax` | maximum subsurface flow | m``^2`` Δt``{-1}``  | - |
| `to_river` | part of subsurface flow that flows to the river | m``^3`` Δt``{-1}``  | - |
| `wb_pit` | boolean location (0 or 1) of a waterbody (wb, reservoir or lake) | -  | - |