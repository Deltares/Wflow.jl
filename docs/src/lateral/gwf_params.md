# Groundwater flow

## Confined aquifer
The Table below shows the parameters (fields) of struct `ConfinedAquifer`, including a
description of these parameters, the unit, and default value if applicable. Struct
`ConfinedAquifer` is not (yet) part of a Model in Wflow.

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

## Unconfined aquifer
The Table below shows the parameters (fields) of struct `UnconfinedAquifer`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static input data (netCDF), and
can be listed in the TOML configuration file under `[lateral.subsurface]`, to map the
internal model parameter to the external netCDF variable. For some input parameters the
parameter listed under `[lateral.subsurface]` is not equal to the internal model parameter,
these are listed in the Table below between parentheses after the internal model parameter.

|  parameter  | description  	        | unit  | default | 
|:--------------- | ------------------| ----- | -------|
| **`k`** (`conductivity`) | horizontal conductivity  | m d``^{-1}``s | - |
| **`specific_yield`**    | specific yield  | m m``^{-1}`` | - |
| `top`     | top groundwater layers  | m | - |
| `bottom`     | bottom groundwater layers  | m | - |
| `area`          | cell area    | m``^2`` | - |
| `head`          | groundwater head     | m | - |
| `conductance`          | conductance    | m``^2`` d``^{-1}`` | - |

## Constant Head
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

## Boundary conditions

### River
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

### Drainage
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
| `index` (`drain`) |  drain cell index | - | - | 
| `flux`          | exchange flux (drains to aquifer)  | m``^3`` day``^{-1}`` | - |

### Recharge
The Table below shows the parameters (fields) of struct `Recharge`, including a
description of these parameters, the unit, and default value if applicable.

|  parameter  | description  	        | unit  | default |
|:--------------- | ------------------| ----- |
| `rate`          | recharge rate  | m``^3`` day``^{-1}`` | - | - |
| `index`         |  recharge cell index | - | - | - |
| `flux`          | recharge flux  | m``^3`` day``^{-1}`` | 

### Head boundary
The Table below shows the parameters (fields) of struct `HeadBoundary`, including a
description of these parameters, the unit, and default value if applicable.

|  parameter  | description  	        | unit  | default |
|:--------------- | ------------------| ----- | ---- |
| `head`          | head  | m | - |
| `conductance`   |  conductance of the head boundary  | m``^2`` day``^{-1}`` | - |
| `index`         |  head boundary cell index | - | - |
| `flux`          |  conductance of the head boundary  | m``^3`` day``^{-1}`` | - |


### Well boundary
The Table below shows the parameters (fields) of struct `Well`, including a
description of these parameters, the unit, and default value if applicable.

|  input parameter  | description  	        | unit  | default |
|:--------------- | ------------------| ----- | ---- |
| `volumetric_rate` | volumetric well rate  | m``^3`` d``^{-1}`` | - |
| `index` | well index  | - | - |
| `flux`          |  actual well flux  | m``^3`` day``^{-1}`` | - |