## Overland flow
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

## River flow
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