## SBM
The Table below shows the parameters (fields) of struct `SBM`, including a description of
these parameters, the unit, and default value if applicable. The parameters in bold
represent model parameters that can be set through static and forcing input data (netCDF),
and can be listed in the TOML configuration file under `[input.vertical]`, to map the
internal model parameter to the external netCDF variable. For some input parameters the
parameter listed under `[input.vertical]` is not equal to the internal model parameter,
these are listed in the Table below between parentheses after the internal model parameter.
For example, internal model parameter `sl` is mapped as follows in the TOML file to the
external netCDF variable `Sl`:

```toml
[input.vertical]
specific_leaf = "Sl"
```

|  parameter | description    | unit | default |
|:---------------| --------------- | ---------------------- | ----- |
| **`cfmax`**  | degree-day factor | mm ᵒC``^{-1}`` Δt``^{-1}`` | 3.75653 mm ᵒC``^{-1}`` day``^{-1}``  |
| **`tt`**  | threshold temperature for snowfall| ᵒC | 0.0  |
| **`tti`**  | threshold temperature interval length | ᵒC | 1.0  |
| **`ttm`**  | threshold temperature for snowmelt  | ᵒC | 0.0  |
| **`whc`**  | water holding capacity as fraction of current snow pack  | - | 0.1  |
| **`w_soil`**  | soil temperature smooth factor  | - | 0.1125  |
| **`cf_soil`**  | controls soil infiltration reduction factor when soil is frozen  | - | 0.038  |
| **`g_tt`**  | threshold temperature for snowfall above glacier  | ᵒC| 0.0  |
| **`g_cfmax`**  | Degree-day factor for glacier  | mm ᵒC``^{-1}`` Δt``^{-1}``| 3.0 mm ᵒC``^{-1}`` day``^{-1}`` |
| **`g_sifrac`**  | fraction of the snowpack on top of the glacier converted into ice  | - | 0.001  |
| **`glacierfrac`**  | fraction covered by a glacier | - | 0.0  |
| **`glacierstore`**  | water within the glacier  | mm | 5500.0  |
| **`θₛ`** (`theta_s`) | saturated water content (porosity) | - | 0.6  |
| **`θᵣ`** (`theta_r`) | residual water content | - | 0.01  |
| **`kv₀`** (`kv_0`) | Vertical hydraulic conductivity at soil surface | mm Δt``^{-1}`` | 3000.0 mm day``^{-1}``|
| **`f`** | scaling parameter (controls exponential decline of kv₀) | mm``^{-1}`` | 0.001  |
| **`hb`** | air entry pressure of soil (Brooks-Corey) | cm | 10.0  |
| **`soilthickness`** | soil thickness | mm | 2000.0  |
| **`infiltcappath`** | infiltration capacity of the compacted areas  | mm Δt``^{-1}`` | 10.0 mm day``^{-1}`` |
| **`infiltcapsoil`** | soil infiltration capacity | mm Δt``^{-1}`` | 100.0 mm day``^{-1}``|
| **`maxleakage`** | maximum leakage from saturated zone | mm Δt``^{-1}`` | 0.0  mm day``^{-1}``|
| **`c`** | Brooks-Corey power coefficient for each soil layer  | - | 10.0  |
| **`kvfrac`** | muliplication factor applied to kv_z (vertical flow)  | - | 1.0  |
| **`waterfrac`** | fraction of open water (excluding rivers) | - | 0.0  |
| **`pathfrac`** | fraction of compacted area | - | 0.01  |
| **`rootingdepth`** | rooting depth  | mm | 750.0  |
| **`rootdistpar`** | controls how roots are linked to water table | - | -500.0  |
| **`cap_hmax`** | controlling capillary rise | mm | 2000.0  |
| **`et_reftopot`** | multiplication factor to correct reference evaporation | - | 1.0  |
| **`sl`** (`specific_leaf`) | specific leaf storage  | mm | - |
| **`swood`** (`storage_wood`) | storage woody part of vegetation | mm | - |
| **`kext`** | extinction coefficient (to calculate canopy gap fraction) | - | - |
| **`cmax`** | maximum canopy storage | mm | 1.0 |
| **`e_r`** (`eoverr`) | Gash interception model parameter | - | 0.1 |
| **`canopygapfraction`** | canopy gap fraction | - | 0.1 | - |
| `Δt`             | model time step     | s | - |
| `maxlayers`      | maximum number of soil layers     | - | - |
| `n`      |  number of grid cells    | - | - |
| `nlayers`      |  number of soil layers    | - | - |
| `n_unsatlayers`      |  number of unsaturated soil layers    | - | - |
| `riverfrac`      |  fraction of river    | - | - |
| `act_thickl`      |  thickness of soil layers    | mm | - |
| `sumlayers`      |  cumulative sum of soil layers thickness, starting at soil surface | mm | - |
| `stemflow`|  stemflow | mm Δt``^{-1}`` | - |
| `throughfall`|  throughfall | mm Δt``^{-1}`` | - |
| `ustorelayerdepth`|  amount of water in the unsaturated store, per layer | mm | - |
| `satwaterdepth`|  saturated store | mm | - |
| `zi`|  pseudo-water table depth (top of the saturated zone) | mm | - |
| `soilwatercapacity`|  soilwater capacity | mm | - |
| `canopystorage`|  canopy storage | mm | - |
| `canopygapfraction` | canopygapfraction | - | - |
|**`precipitation`** | precipitation | mm Δt``^{-1}``|  - |
| **`temperature`** | temperature | ᵒC | - |
| **`potential_evaporation`** | potential evaporation | mm Δt``^{-1}`` | - |
| `pottrans_soil` | interception subtracted from potential evaporation) | mm Δt``^{-1}`` | - |
| `transpiration` | transpiration | mm Δt``^{-1}`` | - |
| `ae_ustore` | actual evaporation from unsaturated store | mm Δt``^{-1}`` | - |
| `interception` | interception | mm Δt``^{-1}`` | - |
| `soilevap` | soil evaporation from unsaturated store | mm Δt``^{-1}`` | - |
| `soilevapsat` | soil evaporation from saturated store | mm Δt``^{-1}`` | - |
| `actcapflux` | actual capillary rise | mm Δt``^{-1}`` | - |
| `actevapsat` | actual transpiration from saturated store | mm Δt``^{-1}`` | - |
| `actevap` | total actual evapotranspiration  | mm Δt``^{-1}`` | - |
| `runoff_river` | runoff from river based on `riverfrac`  | mm Δt``^{-1}`` | - |
| `runoff_land` | runoff from land based on `waterfrac`  | mm Δt``^{-1}`` | - |
| `ae_openw_l` | actual evaporation from open water (land)  | mm Δt``^{-1}`` | - |
| `ae_openw_r` | actual evaporation from river  | mm Δt``^{-1}`` | - |
| `net_runoff_river` | net runoff from river (`runoff_river` - `ae_openw_r`)  | mm Δt``^{-1}`` | - |
| `avail_forinfilt` | water available for infiltration  | mm Δt``^{-1}`` | - |
| `actinfilt` | actual infiltration into the unsaturated zone  | mm Δt``^{-1}`` | - |
| `actinfiltsoil` | actual infiltration into non-compacted fraction  | mm Δt``^{-1}`` | - |
| `actinfiltpath` | actual infiltration into compacted fraction  | mm Δt``^{-1}`` | - |
| `infiltsoilpath` | infiltration into the unsaturated zone | mm Δt``^{-1}`` | - |
| `infiltexcess` | infiltration excess water | mm Δt``^{-1}`` | - |
| `excesswater` | water that cannot infiltrate due to saturated soil (saturation excess) | mm Δt``^{-1}`` | - |
| `exfiltsatwater` | water exfiltrating during saturation excess conditions | mm Δt``^{-1}`` | - |
| `exfiltustore` | water exfiltrating from unsaturated store because of change in water table | mm Δt``^{-1}`` | - |
| `excesswatersoil` | excess water for non-compacted fraction | mm Δt``^{-1}`` | - |
| `excesswaterpath` | excess water for ompacted fraction | mm Δt``^{-1}`` | - |
| `runoff` | total surface runoff from infiltration and saturation excess  | mm Δt``^{-1}`` | - |
| `vwc` | volumetric water content per soil layer (including θᵣ and saturated zone)  | - | - |
| `vwc_perc` | volumetric water content per soil layer (including θᵣ and saturated zone)  | % | - |
| `rootstore` | root water storage in unsaturated and saturated zone (excluding θᵣ)  | mm| - |
| `vwc_root` | volumetric water content in root zone (including θᵣ and saturated zone) | -| - |
| `vwc_percroot` | volumetric water content in root zone (including θᵣ and saturated zone) | % | - |
| `ustoredepth` | total amount of available water in the usaturated zone | mm | - |
| `transfer` | downward flux from unsaturated to saturated zone | mm Δt``^{-1}`` | - |
| `recharge` | net recharge to saturated zone | mm Δt``^{-1}`` | - |
| `actleakage` | actual leakage from saturated store  | mm Δt``^{-1}`` | - |
| `snow` | snow storage  | mm | - |
| `snowwater` | liquid water content in the snow pack  | mm | - |
| `rainfallplusmelt` | snowmelt + precipitation as rainfall | mm Δt``^{-1}`` | - |
| `glacierstore` | water within the glacier | mm | - |
| `tsoil` | top soil temperature | ᵒC | - |
| **`leaf_area_index`** | leaf area index | m``^2`` m``{-2}`` | - |
| `waterlevel_land` | water level land | mm | - |
| `waterlevel_river` | water level river | mm | - |


## HBV
The Table below shows the parameters (fields) of struct `HBV`, including a description of
these parameters, the unit, and default value if applicable. The parameters in bold
represent model parameters that can be set through static and forcing input data (netCDF),
and can be listed in the TOML configuration file under `[input.vertical]`, to map the
internal model parameter to the external netCDF variable. 

|  parameter | description    | unit | default |
|:---------------| --------------- | ---------------------- | ----- |
| **`cfmax`**  | degree-day factor | mm ᵒC``^{-1}`` Δt``^{-1}`` | 3.75653 mm ᵒC``^{-1}`` day``^{-1}``  |
| **``tt`**  | threshold temperature for snowfall| ᵒC | -1.41934  |
| **`tti`**  | threshold temperature interval length | ᵒC | 1.0  |
| **`ttm`**  | threshold temperature for snowmelt  | ᵒC | -1.41934  |
| **`whc`**  | water holding capacity as fraction of current snow pack  | - | 0.1  |
| **`g_tt`**  | threshold temperature for snowfall above glacier  | ᵒC| 0.0  |
| **`g_cfmax`**  | Degree-day factor for glacier  | mm ᵒC``^{-1}`` Δt``^{-1}``| 3.0 mm ᵒC``^{-1}`` day``^{-1}`` |
| **`g_sifrac`**  | fraction of the snowpack on top of the glacier converted into ice  | - | 0.001  |
| **`glacierfrac`**  | fraction covered by a glacier | - | 0.0  |
| **`glacierstore`**  | water within the glacier  | mm | 5500.0  |
| **`fc`**  | field capacity  | mm | 260.0  |
| **`betaseepage`**  | exponent in soil runoff generation equation  | - | 1.8  |
| **`lp`**  | fraction of field capacity below which actual evaporation=potential evaporation | - | 0.53 |
| **`k4`**  | recession constant baseflow | Δt``^-1`` | 0.02307 day``^{-1}`` |
| **`kquickflow`**  | recession constant upper reservoir | Δt``^-1`` |  0.09880 day``^{-1}`` |
| **`suz`**  | Level over which `k0` is used | mm |  100.0 |
| **`k0`**  | recession constant upper reservoir | Δt``^-1`` | 0.30 day``^{-1}`` |
| **`khq`**  | recession rate at flow `hq`| Δt``^-1`` | 0.09880 day``^{-1}`` |
| **`hq`**  | high flow rate hq for which recession rate of upper reservoir is known | mm Δt``^-1`` | 3.27 mm day``^{-1}`` |
| **`alphanl`**  | measure of non-linearity of upper reservoir | - | 1.1 |
| **`perc`**  | percolation from upper to lower zone |  mm Δt``^-1`` | 0.4 mm day``^{-1}`` |
| **`cfr`**  | refreezing efficiency constant in refreezing of freewater in snow | - | 0.05  |
| **`pcorr`**  | correction factor for precipitation  | - | 1.0  |
| **`rfcf`**  | correction factor for rainfall   | - | 1.0  |
| **`sfcf`**  | correction factor for snowfall | - | 1.0  |
| **`cflux`**  | maximum capillary rise from runoff response routine to soil moisture routine | mm Δt``^-1`` | 2.0 mm day``^{-1}``  |
| **`icf`**  | maximum interception storage (in forested and non-forested areas) | mm | 2.0  |
| **`cevpf`**  | correction factor for potential evaporation | - | 1.0  |
| **`epf`**  | exponent of correction factor for evaporation on days with precipitation | mm``^{-1}`` | 1.0  |
| **`ecorr`**  | evaporation correction | - | 1.0  |
| **`precipitation`**           | precipitation     | mm Δt``^-1`` | - |
| **`temperature`**             | temperature    | ᵒC | - |
| **`potential_evaporation`**            | potential evapotranspiration     | mm Δt``^-1`` | - |
| `potsoilevap`             | potential soil evaporation      | mm Δt``^-1`` | - |
| `soilevap`             | soil evaporation | mm Δt``^-1`` | - |
| `intevap`             | evaporation from interception storage     | mm Δt``^-1`` | - |
| `actevap`             | actual evapotranspiration (intevap + soilevap) | mm Δt``^-1`` | - |
| `interceptionstorage`             | actual interception storage     | mm | - |
| `snowwater`             | available free water in snow    | mm | - |
| `snow`             | snow pack    | mm | - |
| `rainfallplusmelt`    | snow melt + precipitation as rainfall    | mm Δt``^-1`` | - |
| `soilmoisture`        | actual soil moisture | mm | - |
| `directrunoff`        | direct runoff to upper zone    | mm Δt``^-1`` | - |
| `hbv_seepage`         | recharge to upper zone  | mm Δt``^-1``  | - |
| `in_upperzone`         | water inflow into upper zone  | mm Δt``^-1``  | - |
| `upperzonestorage`     | water content of the upper zone  | mm  | - |
| `quickflow`         |  specific runoff (quickflow part)  | mm Δt``^-1``  | - |
| `real_quickflow`         | specific runoff (quickflow), if K upper zone is precalculated  | mm Δt``^-1``  | - |
| `percolation`         | actual percolation to the lower zone | mm Δt``^-1``  | - |
| `capflux`         | capillary rise | mm Δt``^-1``  | - |
| `lowerzonestorage`         | water content of the lower zone | mm  | - |
| `baseflow`         | specific runoff (baseflow part) per cell  | mm Δt``^-1``  | - |
| `runoff`         | total specific runoff per cell  | mm Δt``^-1``  | - |


## Sediment

The Table below shows external parameters that can be set through static input data
(netCDF), and can be listed in the TOML configuration file under `[input.vertical]`. These
external parameters are not part of struct `LandSediment`, but used to calculate parameters
of struct `LandSediment`.

| external parameter | description    | unit | default |
|:---------------| --------------- | ---------------------- | ----- |
| `pclay` | percentage clay | %  | 0.1 |
| `psilt` | percentage silt  | % | 0.1 |
| `resareas` | reservoir coverage | - | - |
| `lakeareas` | lake coverage | - | - |

The Table below shows the parameters (fields) of struct `LandSediment`, including a
description of these parameters, the unit, and default value if applicable. The parameters
in bold represent model parameters that can be set through static and forcing input data
(netCDF), and can be listed in the TOML configuration file under `[input.vertical]`, to map
the internal model parameter to the external netCDF variable. For some input parameters the
parameter listed under `[input.vertical]` is not equal to the internal model parameter,
these are listed in the Table below between parentheses after the internal model parameter.
For example, internal model parameter `sl` is mapped as follows in the TOML file to the
external netCDF variable `Sl`:

```toml
[input.vertical]
specific_leaf = "Sl"
```

|  parameter | description    | unit | default |
|:---------------| --------------- | ---------------------- | ----- |
| **`canopyheight`**  | canopy height | m | 3.0  |
| **`erosk`**  | coefficient for EUROSEM rainfall erosion | - | 0.6  |
| **`erosspl`**  | exponent for EUROSEM rainfall erosion | - | 2.0  |
| **`erosov`**  | coefficient for ANSWERS overland flow erosion | - | 0.9  |
| **`pathfrac`**  | fraction of impervious area per grid cell | - | 0.01  |
| **`slope`**  | land slope | - | 0.01  |
| **`usleC`**  | USLE crop management factor | - | 0.01  |
| **`usleK`**  | USLE soil erodibility factor | - | 0.1 |
| **`sl`** (`specific_leaf`) | specific leaf storage  | mm | - |
| **`swood`** (`storage_wood`) | storage woody part of vegetation | mm | - |
| **`kext`** | extinction coefficient (to calculate canopy gap fraction) | - | - |
| **`cmax`** | maximum canopy storage | mm | 1.0 |
| **`canopygapfraction`** | canopy gap fraction | - | 0.1 |
| **`dmclay`** | median diameter particle size class clay | µm | 2.0 |
| **`dmsilt`** | median diameter particle size class silt | µm| 10.0 |
| **`dmsand`** | median diameter particle size class sand | µm | 200.0 |
| **`dmsagg`** | median diameter particle size class small aggregates | µm | 30.0 |
| **`dmlagg`** | median diameter particle size class large aggregates | µm | 500.0 |
| **`rhos`** (`rhosed`) | density of sediment | kg m``^{-3}1`` | 2650.0 |
| `n`             | number of cells     | - | - |
| `yl`            | length of cells in y direction | m | - |
| `xl`            | length of cells in x direction     | m | - |
| `riverfrac`     | fraction of river   | - | - |
| `wbcover`     | waterbody coverage   | - | - |
| **`h_land`**     | depth of overland flow   | m | - |
| **`interception`**     | canopy interception   | mm Δt``^{-1}`` | - |
| **`precipitation`**     | precipitation  | mm Δt``^{-1}`` | - |
| **`q_land`**     | overland flow  | m``^3`` s``^{-1}`` | - |
| `sedspl`     | sediment eroded by rainfall  | ton Δt``^{-1}`` | - |
| `sedov`     | sediment eroded by overland flow  | ton Δt``^{-1}`` | - |
| `soilloss`     | total eroded soil  | ton Δt``^{-1}`` | - |
| `erosclay`     | eroded soil for particle class clay  | ton Δt``^{-1}`` | - |
| `erossilt`     | eroded soil for particle class silt  | ton Δt``^{-1}`` | - |
| `erossand`     | eroded soil for particle class sand  | ton Δt``^{-1}`` | - |
| `erossagg`     | eroded soil for particle class small aggregates | ton Δt``^{-1}`` | - |
| `eroslagg`     | eroded soil for particle class large aggregates | ton Δt``^{-1}`` | - |
| **`leaf_area_index`** | leaf area index | m``^2`` m``^{-2}`` | - |
| `dl` | drain length | m | - |
| `dw` | flow width | m | - |
| `cGovers` | Govers transport capacity coefficient | - | - |
| `nGovers` | Govers transport capacity coefficient | - | - |
| `D50` | median particle diameter of the topsoil | mm | - |
| `fclay` | fraction of particle class clay | - | - |
| `fsilt` | fraction of particle class silt | - | - |
| `fsand` | fraction of particle class sand | - | - |
| `fsagg` | fraction of particle class small aggregates | - | - |
| `flagg` | fraction of particle class large aggregates  | - | - |
| `rivcell` | river cells  | - | - |
| `TCsed` | total transport capacity of overland flow  | ton Δt``^{-1}`` | - |
| `TCclay` | transport capacity of overland flow for particle class clay | ton Δt``^{-1}``| - |
| `TCsilt` | transport capacity of overland flow for particle class silt  | ton Δt``^{-1}`` | - |
| `TCsand` | transport capacity of overland flow for particle class sand  | ton Δt``^{-1}`` | - |
| `TCsagg` | transport capacity of overland flow for particle class small aggregates  | ton Δt``^{-1}`` | - |
| `TClagg` | transport capacity of overland flow for particle class large aggregates  | ton Δt``^{-1}`` | - |