# [Parameters lateral concepts](@id params_lat)

## Kinematic wave

### Surface flow
The Table below shows the parameters (fields) of struct `SurfaceFlow`, including a
description of these parameters, the unit, and default value if applicable. `SurfaceFlow` is
used for river and overland flow. The parameters in bold represent model parameters that can
be set through static input data (netCDF), and can be listed in the TOML configuration file
under `[input.lateral.river]` and `[input.lateral.land]`, for river and overland flow
respectively, to map the internal model parameter to the external netCDF variable. For river
flow three additionally parameters can be set, `dl` (river length), `width` (river width)
and bankfull river depth `bankfull_depth` as follows, through the TOML file:

```toml
[input.lateral.river]
length = "wflow_riverlength"
width = "wflow_riverwidth"
bankfull_depth = "wflow_riverdepth"
```
[^1]: default value for Manning's roughness `n`, river = 0.036; land = 0.072
[^2]: only applicable for river domain

|  parameter  | description  	  | unit  | default |
|:--------------- | ------------------| ----- | -------- |
| `β`             |  constant in Manning's equation | - | - |
| `dl`            |  length      | m     |  -      |
| **`n`**             | Manning's roughness | s m``^{-\frac{1}{3}}``| [^1] |
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
| `bankfull_depth`   | bankfull river depth  | m | 1.0 | 
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

### [Reservoirs](@id reservoir_params)
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

### [Lakes](@id lake_params)
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

### Lateral subsurface flow
The Table below shows the parameters (fields) of struct `LateralSSF`, including a
description of these parameters, the unit, and default value if applicable. 

| parameter | description  	        | unit | default |
|:---------------| --------------- | ---------------------- | ----- |    
| `kh₀`          | horizontal hydraulic conductivity at soil surface  | m d``^{-1}`` |  - |
| `f` | a scaling parameter (controls exponential decline of kh₀) | m``^{-1}``  | - |
| `soilthickness` | soil thickness | m  | - |
| `θₛ` | saturated water content (porosity) | -  | - |
| `θᵣ` | residual water content  | -  | - |
| `Δt` | model time step | d  | - |
| `βₗ` | slope | -  | - |
| `dl` | drain length | m | - |
| `dw` | drain width | m | - |
| `zi` | pseudo-water table depth (top of the saturated zone) | m | - |
| `exfiltwater` | exfiltration (groundwater above surface level, saturated excess conditions) | m Δt⁻¹ | - |
| `recharge` | net recharge to saturated store  | m Δt⁻¹ | - |
| `ssf` | subsurface flow  | m``^3`` d``{-1}``  | - |
| `ssfin` | inflow from upstream cells | m``^3`` d``{-1}``  | - |
| `ssfmax` | maximum subsurface flow | m``^2`` d``{-1}``  | - |
| `to_river` | part of subsurface flow that flows to the river | m``^3`` d``{-1}``  | - |
| `wb_pit` | boolean location (0 or 1) of a waterbody (wb, reservoir or lake) | -  | - |

## Local inertial

### [River flow](@id local-inertial_river_params)
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

### [Overland flow](@id local-inertial_land_params)
The Table below shows the parameters (fields) of struct `ShallowWaterLand`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[input.lateral.land]`, to map the
internal model parameter to the external netCDF variable.

[^3]: for river cells the reference level is the river bed elevation `zb`

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
| `h` [^3] |  water depth of cell | m | - |
| **`z`**  |  elevation of cell | m | - |
| `froude_limit`  |  if true a check is performed if froude number > 1.0 (algorithm is modified)| - | - |
| `rivercells`  |  river cells| - | - |
| `h_av` [^3]  | average water depth| m | - |

## Groundwater flow

### Confined aquifer
The Table below shows the parameters (fields) of struct `ConfinedAquifer`, including a
description of these parameters, the unit, and default value if applicable. Struct
`ConfinedAquifer` is not (yet) part of a Wflow Model.

|  parameter  | description  	      | unit  | default | 
|:--------------- | ------------------| ----- | -------|
| `k` | horizontal conductivity  | m d``^{-1}``s | - |
| `storativity`     | storativity  | m m``^{-1}`` | - |
| `specific_storage` | specific storage | m``^{-1}`` | - }
| `top`     | top groundwater layers  | m | - |
| `bottom`     | bottom groundwater layers  | m | - |
| `area`          | cell area    | m``^2`` | - |
| `head`          | groundwater head     | m | - |
| `conductance`          | conductance    | m``^2`` d``^{-1}`` | - |

### Unconfined aquifer
The Table below shows the parameters (fields) of struct `UnconfinedAquifer`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[lateral.subsurface]`, to map the
internal model parameter to the external netCDF variable. For some input parameters the
parameter listed under `[lateral.subsurface]` is not equal to the internal model parameter,
these are listed in the Table below between parentheses after the internal model parameter.
The `top` parameter is provided by the external parameter `altitude` as part of the static
input data and set as follows through the TOML file:

```toml
[input]
# these are not directly part of the model
altitude = "wflow_dem"
```

|  parameter  | description  	        | unit  | default | 
|:--------------- | ------------------| ----- | -------|
| **`k`** (`conductivity`) | horizontal conductivity  | m d``^{-1}``s | - |
| **`specific_yield`**    | specific yield  | m m``^{-1}`` | - |
| **`top`** (`altitude`)  | top groundwater layer  | m | - |
| `bottom`     | bottom groundwater layer  | m | - |
| `area`          | cell area    | m``^2`` | - |
| `head`          | groundwater head     | m | - |
| `conductance`          | conductance    | m``^2`` d``^{-1}`` | - |

### Constant Head
The Table below shows the parameters (fields) of struct `ConstantHead`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[lateral.subsurface]`, to map the
internal model parameter to the external netCDF variable. For some input parameters the
parameter listed under `[lateral.subsurface]` is not equal to the internal model parameter,
these are listed in the Table below between parentheses after the internal model parameter.

|  parameter  | description  	       | unit  | default |
|:--------------- | ------------------| ----- | --------- |
| **`head`** (`constant_head`)          | groundwater head     | m | - |
| `index`          | constand head cell index | - | - |

### Boundary conditions

#### [River](@id gwf_river_params)
The Table below shows the parameters (fields) of struct `River`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[lateral.subsurface]`, to map the
internal model parameter to the external netCDF variable. For some input parameters the
parameter listed under `[lateral.subsurface]` is not equal to the internal model parameter,
these are listed in the Table below between parentheses after the internal model parameter.

|  parameter  | description  	        | unit  | default | 
|:--------------- | ------------------| ----- | -------|
| `stage`    | river stage  | m | - |
| **`infiltration_conductance`**     | river bed infiltration conductance  | m``^2`` day``^{-1}`` m``^2`` day``^{-1}``| - |
| **`exfiltration_conductance`**     | river bed exfiltration conductance   | m``^2`` day``^{-1}`` | - |
| **`bottom`** (`river_bottom`)     | river bottom elevation  | m | - |
| `index`         |  river cell index | - | - |
| `flux`          | exchange flux (river to aquifer)  | m``^3`` d``^{-1}`` | - |

#### [Drainage](@id gwf_drainage_params)
The Table below shows the parameters (fields) of struct `Drainage`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[lateral.subsurface]`, to map the
internal model parameter to the external netCDF variable. For some input parameters the
parameter listed under `[lateral.subsurface]` is not equal to the internal model parameter,
these are listed in the Table below between parentheses after the internal model parameter.

|  parameter  | description  	        | unit  | default | 
|:--------------- | ------------------| ----- | -------|
| **`elevation`** (`drain_elevation`)   | drain elevation  | m | - |
| **`conductance`** (`drain_conductance`)  | drain conductance  | m``^2`` day``^{-1}`` | - |
| **`index`** (`drain`) |  drain cell index | - | - | 
| `flux`          | exchange flux (drains to aquifer)  | m``^3`` day``^{-1}`` | - |

#### [Recharge](@id gwf_recharge_params)
The Table below shows the parameters (fields) of struct `Recharge`, including a
description of these parameters, the unit, and default value if applicable.

|  parameter  | description  	        | unit  | default |
|:--------------- | ------------------| ----- | ---- |
| `rate`          | recharge rate  | m``^3`` day``^{-1}`` | - |
| `index`         |  recharge cell index | - | - |
| `flux`          | recharge flux  | m``^3`` day``^{-1}`` | - |

#### [Head boundary](@id gwf_headboundary_params)
The Table below shows the parameters (fields) of struct `HeadBoundary`, including a
description of these parameters, the unit, and default value if applicable.

|  parameter  | description  	        | unit  | default |
|:--------------- | ------------------| ----- | ---- |
| `head`          | head  | m | - |
| `conductance`   |  conductance of the head boundary  | m``^2`` day``^{-1}`` | - |
| `index`         |  head boundary cell index | - | - |
| `flux`          |  conductance of the head boundary  | m``^3`` day``^{-1}`` | - |


#### [Well boundary](@id well_boundary_params)
The Table below shows the parameters (fields) of struct `Well`, including a description of
these parameters, the unit, and default value if applicable.

|  input parameter  | description  	        | unit  | default |
|:--------------- | ------------------| ----- | ---- |
| `volumetric_rate` | volumetric well rate  | m``^3`` d``^{-1}`` | - |
| `index` | well index  | - | - |
| `flux`          |  actual well flux  | m``^3`` day``^{-1}`` | - |

## Sediment

### Overland flow
The Table below shows the parameters (fields) of struct `OverlandFlowSediment`, including a
description of these parameters, the unit, and default value if applicable. 

|  parameter | description  	        | unit | default |
|:---------------| --------------- | ---------------------- | ------- |
| `n`             | number of cells     | - | - |
| `rivcell`       | river cells | - | - |
| `soilloss`       | total eroded soil | ton Δt``^{-1}`` | - |
| `erosclay`       | eroded soil for particle class clay | ton Δt``^{-1}`` | - |
| `erossilt`       | eroded soil for particle class silt | ton Δt``^{-1}`` | - |
| `erossand`       | eroded soil for particle class sand | ton Δt``^{-1}`` | - |
| `erossagg`       | eroded soil for particle class small aggregates | ton Δt``^{-1}`` | - |
| `eroslagg`       | eroded soil for particle class large aggregates | ton Δt``^{-1}`` | - |
| `TCsed`        | total transport capacity of overland flow | ton Δt``^{-1}`` | - |
| `TCclay`       | transport capacity of overland flow for particle class clay | ton Δt``^{-1}`` | - |
| `TCsilt`       | transport capacity of overland flow for particle class silt | ton Δt``^{-1}`` | - |
| `TCsand`       | transport capacity of overland flow for particle class sand | ton Δt``^{-1}`` | - |
| `TCsagg`       | transport capacity of overland flow for particle class small aggregates | ton Δt``^{-1}`` | - |
| `TClagg`       | transport capacity of overland flow for particle class large aggregates | ton Δt``^{-1}`` | - |
| `inlandsed`       | sediment reaching the river with overland flow | ton Δt``^{-1}`` | - |
| `inlandclay`       | sediment with particle class clay reaching the river with overland flow | ton Δt``^{-1}`` | - |
| `inlandsilt`       | sediment with particle class silt reaching the river with overland flow | ton Δt``^{-1}`` | - |
| `inlandsand`       | sediment with particle class sand reaching the river with overland flow | ton Δt``^{-1}`` | - |
| `inlandsagg`       | sediment with particle class small aggregates reaching the river with overland flow | ton Δt``^{-1}`` | - |
| `inlandlagg`       | sediment with particle class large aggregates reaching the river with overland flow | ton Δt``^{-1}`` | - |

### River flow
The Table below shows external parameters that can be set through static input data
(netCDF), and can be listed in the TOML configuration file under `[input.lateral.river]`.
These external parameters are not part of struct `RiverSediment`, but used to calculate
parameters of struct `RiverSediment`.

| external parameter | description    | unit | default |
|:---------------| --------------- | ---------------------- | ----- |
| `reslocs`  | reservoir location (outlet) | - | -  |
| `resareas`  | reservoir coverage | - | -  |
| `resarea`  | reservoir area | - | m``^2``  |
| `restrapeff`  | reservoir trapping efficiency coefficient | - | -  |
| `lakelocs`  | lake location (outlet) | - | -  |
| `lakeareas`  | lake coverage | - | -  |
| `lakearea`  | lake area | - | m``^2``  | - |


The Table below shows the parameters (fields) of struct `RiverSediment`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static and forcing input data
(netCDF), and can be listed in the TOML configuration file under `[input.lateral.river]`, to
map the internal model parameter to the external netCDF variable. For some input parameters
the parameter listed under `[input.lateral.river]` is not equal to the internal model
parameter, these are listed in the Table below between parentheses after the internal model
parameter. For example, internal model parameter `sl` is mapped as follows in the TOML file
to the external netCDF variable `RiverSlope`:

```toml
[input.vertical]
slope = "RiverSlope"
```

| parameter | description    | unit | default |
|:---------------| --------------- | ---------------------- | ----- |
| **`dl`** (`length`)            | river length    | m | - |
| **`width`**           | river width    | m | - |
| **`sl`** (`slope`)  | river slope | - | -  |
| **`rhos`** (`rhosed`)  | density of sediment | kg m``^{-3}1`` | 2650.0 |
| **`dmclay`** | median diameter particle size class clay | mm | 2.0 |
| **`dmsilt`** | median diameter particle size class silt | mm| 10.0 |
| **`dmsand`** | median diameter particle size class sand | mm | 200.0 |
| **`dmsagg`** | median diameter particle size class small aggregates | mm | 30.0 |
| **`dmlagg`** | median diameter particle size class large aggregates | mm | 500.0 |
| **`dmgrav`** | median diameter particle size class gravel | mm | 2000.0 |
| **`fclayriv`** | fraction of particle class clay | - | - |
| **`fsiltriv`** | fraction of particle class silt | - | - |
| **`fsandriv`** | fraction of particle class sand | - | - |
| **`fsaggriv`** | fraction of particle class small aggregates | - | - |
| **`flaggriv`** | fraction of particle class large aggregates  | - | - |
| **`fgravriv`** | fraction of particle class gravel  | - | - |
| **`d50`** (`d50riv`) | river sediment median diameter  | mm | - |
| **`d50engelund`**| river mean diameter  | mm | - |
| **`cbagnold`**| Bagnold c coefficient  | - | - |
| **`ebagnold`**| Bagnold exponent | - | - |
| `n`             | number of cells     | - |  - |
| `Δt`             | model time step     | s |  - |
| `ak`             | Kodatie coefficient `a`    | - |  - |
| `bk`             | Kodatie coefficient `b`    | - |  - |
| `ck`             | Kodatie coefficient `c`    | - |  - |
| `dk`             | Kodatie coefficient `d`    | - |  - |
| `kdbank`             | bank erodibilty    | m``^3`` N``^{-1}`` s``^{-1}`` | - |
| `kdbed`             | bed erodibility    | m``^3`` N``^{-1}`` s``^{-1}`` | - |
| `TCrbank`             | critical bed bank shear stress    | m``^3`` N``^{-2}`` | - |
| `TCrbed`             | critical bed shear stress    | m``^3`` N``^{-2}`` | - |
| **`h_riv`**             | river water level  | m| - |
| **`q_riv`**            | river discharge  | m``^3`` s``^{-1}`` | - |
| `inlandclay`             | sediment input with particle class clay from land erosion  | t Δt``^{-1}`` | - |
| `inlandsilt`             | sediment input with particle class silt from land erosion  | t Δt``^{-1}`` | - |
| `inlandsand`             | sediment input with particle class sand from land erosion  | t Δt``^{-1}`` | - |
| `inlandsagg`             | sediment input with particle class small aggregates from land erosion  | t Δt``^{-1}`` | - |
| `inlandlagg`             | sediment input with particle class large aggregates from land erosion  | t Δt``^{-1}`` | - |
| `inlandsed`             | sediment input from land erosion  | t Δt``^{-1}`` | - |
| `sedload`             | sediment left in the cell  | t | - |
| `clayload`             | sediment with particle class clay left in the cell  | t | - |
| `siltload`             | sediment with particle class silt left in the cell  | t  | - |
| `sandload`             | sediment with particle class sand left in the cell  | t  | - |
| `saggload`             | sediment with particle class small aggregates left in the cell  | t  | - |
| `laggload`             | sediment with particle class large aggregates in the cell  | t  | - |
| `gravload`             | sediment with particle class gravel left in the cell  | t  | - |
| `sedstore`              | sediment stored on the river bed after deposition  | t Δt``^{-1}``| - |
| `claystore`             | sediment with particle class clay stored on the river bed after deposition  | t Δt``^{-1}`` | - |
| `siltstore`             | sediment with particle class silt stored on the river bed after deposition | t Δt``^{-1}`` | - |
| `sandstore`             | sediment with particle class sand stored on the river bed after deposition | t Δt``^{-1}`` | - |
| `saggstore`             | sediment with particle class small aggregates stored on the river bed after deposition  | t Δt``^{-1}`` | - |
| `laggstore`             | sediment with particle class large aggregates stored on the river bed after deposition  | t Δt``^{-1}`` | - |
| `gravstore`             | sediment with particle class gravel stored on the river bed after deposition  | t Δt``^{-1}``| - |
| `outsed`              | sediment flux   | t Δt``^{-1}``| - |
| `outclay`              | sediment with particle class clay flux  | t Δt``^{-1}``| - |
| `outsilt`              | sediment with particle class silt  | t Δt``^{-1}``| - |
| `outsand`              | sediment with particle class sand  | t Δt``^{-1}``| - |
| `outsagg`              | sediment with particle class small aggregates  | t Δt``^{-1}``| - |
| `outlagg`              | sediment with particle class large aggregates  | t Δt``^{-1}``| - |
| `outgrav`              | sediment with particle class gravel  | t Δt``^{-1}``| - |
| `Sedconc`              | sediment concentration  | kg m``^{-3}``| - |
| `SSconc`               | sediment concentration   | kg m``^{-3}``| - |
| `Bedconc`              | sediment concentration  | kg m``^{-3}``| - |
| `maxsed`              | river transport capacity | t Δt``^{-1}``| - |
| `erodsed`              | total eroded sediment | t Δt``^{-1}``| - |
| `erodsedbank`          | eroded bank sediment | t Δt``^{-1}``| - |
| `erodsedbed`           | eroded bed sediment  | t Δt``^{-1}``| - |
| `depsed`           | deposited sediment  | t Δt``^{-1}``| - |
| `insed`           | sediment input flux  | t Δt``^{-1}``| - |
| `wbcover`           | waterbody coverage  | - | - |
| `wblocs`           | waterbody locations  | - | - |
| `wbarea`           | waterbody area  | m``^2`` | - |
| `wbtrap`           | waterbody trapping efficiency coefficient  | - |  - |
