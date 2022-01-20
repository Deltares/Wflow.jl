# Config and TOML

[TOML](https://github.com/toml-lang/toml) is used as the configuration format for models
available in wflow. File paths included in the configuration TOML file are relative to the
TOML file location, or to `dir_input` and `dir_output` if they are given.

## General time info
Time information is optional. When left out, each time step in the forcing NetCDF will be
calculated. If you wish to calculate a subset of this time range, or a different timestep,
you can specify a `starttime`, `endtime` and `timestepsecs` yourself. The `time_units`
optional information is used by the `writer` of the model, for model output in netCDF
format. The `calendar` option allows you to calculate in one of the different [CF
conventions calendars](http://cfconventions.org/cf-conventions/cf-conventions.html#calendar)
provided by the [CFTime.jl package](https://juliageo.org/CFTime.jl/latest/), such as
`"360_day"`. This is useful if you want to calculate climate scenarios which are sometimes
provided in these alternative calendars.

```toml
calendar = "standard"                           # optional, this is default value
endtime = 2000-02-01T00:00:00                   # optional, default from forcing NetCDF
starttime = 2000-01-01T00:00:00                 # optional, default from forcing NetCDF
time_units = "days since 1900-01-01 00:00:00"   # optional, this is default value
timestepsecs = 86400                            # optional, default from forcing NetCDF
dir_input = "data/input"                        # optional, default is the path of the TOML
dir_output = "data/output"                      # optional, default is the path of the TOML
```

## Logging
Wflow emits logging messages at various levels such as debug, info, and error. These get
sent to both the terminal as well as a log file. If you want to debug an issue it can be
helpful to set `loglevel = "debug"` in the TOML. To avoid flooding the screen, debug
messages are only sent to the log file. The following settings will affect the logging:

```toml
silent = false          # optional, default is "false"
loglevel = "debug"      # optional, default is "info"
path_log = "log.txt"    # optional, default is "log.txt"
fews_run = false        # optional, default value is false
```

`silent` avoids logging to the terminal, and only writes the log file. `loglevel` controls
which levels are filtered out, so the default setting `"info"` does not show any debug level
messages. Note that for finer control, you can also pass an integer log level, see Julia's
[Logging](https://docs.julialang.org/en/v1/stdlib/Logging/#Log-event-structure)
documentation. `path_log` sets the desired output path for the log file. For information
regarding `fews_run`, see [Run from Delft-FEWS](@ref).

## Model section
Model specific settings can be included in the model section of the TOML file.

```toml
[model]
type = "sbm"                        # one of ("sbm", "sbm_gwf, "hbv")
masswasting = true                  # include lateral snow transport in the model, default is false
snow = true                         # include snow modelling, default is false
reinit = true                       # cold (reinit = true) or warm state (reinit = false), default is true
reservoirs = true                   # include reservoir modelling, default is false
kin_wave_iteration = true           # enable kinematic wave iterations in the model, default is false
thicknesslayers = [100, 300, 800]   # specific SBM setting: for each soil layer a thickness [mm] is specified
min_streamorder = 3                 # minimum stream order to delineate subbasins, default is 4 (for multi-threading computing purposes)
```

## State section
The `state` section in the TOML file provides information about the location of the state
file if the model is initialized with a warm state (`path_input`) and to what file the
states are written at the end of the model run (`path_output`). A mapping between external
state names and internal model states is required. This information is specified for each
model component, the `vertical` model and `lateral` model components. In the example below
the `vertical` component represents the SBM concept, and for the `lateral` components there
is a `river` (including `reservoir`), `land` and `subsurface` domain. The internal model
states are listed on the left side, and the external state names are listed on the right
side. Note that `path_input` is only required when `reinit` is set to false. `path_output`
is optional, an output state file is only written when it is defined. If neither is set,
the entire `state` section can be left out.

```toml
[state]
path_input = "instates-moselle.nc" 
path_output = "outstates-moselle.nc"

[state.vertical]
satwaterdepth = "satwaterdepth"
snow = "snow"
tsoil = "tsoil"
ustorelayerdepth = "ustorelayerdepth"
snowwater = "snowwater"
canopystorage = "canopystorage"

[state.lateral.river]
q = "q_river"
h = "h_river"
h_av = "h_av_river"

[state.lateral.river.reservoir]
volume = "volume_reservoir"

[state.lateral.subsurface]
ssf = "ssf"

[state.lateral.land]
q = "q_land"
h = "h_land"
h_av = "h_av_land"
```
## Input section
The `input` section of the TOML file contains information about the input forcing and model
parameters files (netCDF format). Forcing is applied to the vertical component of the model,
and needs to be mapped to the external netCDF variable name. `forcing` lists the internal
model forcing parameters, and these are mapped to the external netCDF variables listed under
the section `[input.vertical]`. It is possible to provide cyclic parameters to the model. In
the example below this is done for the internal `vertical.leaf_area_index` model parameter,
that is linked to the external netCDF variable "LAI" variable. 

As for the states, mapping of external model parameters (provided in the example below by
the file staticmaps-moselle.nc) is done per model component. If a model parameter is not
mapped a default value will be used if available.

```toml
[input]
# use "forcing-year-*.nc" if forcing files are split in time
path_forcing = "forcing-moselle.nc"
path_static = "staticmaps-moselle.nc"

# these are not directly part of the model
gauges = "wflow_gauges"
ldd = "wflow_ldd"
river_location = "wflow_river"
subcatchment = "wflow_subcatch"

# specify the internal IDs of the parameters which vary over time
# the external name mapping needs to be below together with the other mappings
forcing = [
"vertical.precipitation",
"vertical.temperature",
"vertical.potential_evaporation",
]

cyclic = ["vertical.leaf_area_index"]

[input.vertical]
altitude = "wflow_dem" 
c = "c" 
cf_soil = "cf_soil" 
cfmax = "Cfmax" 
e_r = "EoverR" 
infiltcappath = "InfiltCapPath" 
infiltcapsoil = "InfiltCapSoil" 
kext = "Kext" 
"kv₀" = "KsatVer" 
leaf_area_index = "LAI"
m = "M" 
maxleakage = "MaxLeakage" 
pathfrac = "PathFrac" 
potential_evaporation = "PET" # forcing
precipitation = "P" # forcing
rootdistpar = "rootdistpar" 
rootingdepth = "RootingDepth" 
soilminthickness = "SoilMinThickness" 
soilthickness = "SoilThickness" 
specific_leaf = "Sl" 
storage_wood = "Swood" 
temperature = "TEMP" # forcing
tt = "TT" 
tti = "TTI" 
ttm = "TTM" 
w_soil = "wflow_soil" 
water_holding_capacity = "WHC" 
waterfrac = "WaterFrac" 
"θᵣ" = "thetaR" 
"θₛ" = "thetaS"

[input.lateral.river]
length = "wflow_riverlength"
n = "N_River"
slope = "RiverSlope"
width = "wflow_riverwidth"

[input.lateral.river.reservoir]
area = "ResSimpleArea"
areas = "wflow_reservoirareas"
demand = "ResDemand"
locs = "wflow_reservoirlocs"
maxrelease = "ResMaxRelease"
maxvolume = "ResMaxVolume"
targetfullfrac = "ResTargetFullFrac"
targetminfrac = "ResTargetMinFrac"

[input.lateral.subsurface]
ksathorfrac = "KsatHorFrac"

[input.lateral.land]
n = "N"
slope = "Slope"
```

## Output NetCDF section

### Grid data
This optional section of the TOML file contains the output netCDF file for writing gridded
model output, including a mapping between internal model parameter components and external
netCDF variables.

```toml
[output]
path = "output_moselle.nc"

[output.vertical]
satwaterdepth = "satwaterdepth"
snow = "snow"
tsoil = "tsoil"
ustorelayerdepth = "ustorelayerdepth"
snowwater = "snowwater"
canopystorage = "canopystorage"

[output.lateral.river]
q = "q_river"
h = "h_river"

[output.lateral.river.reservoir]
volume = "volume_reservoir"

[output.lateral.subsurface]
ssf = "ssf"

[output.lateral.land]
q = "q_land"
h = "h_land"
```

### Scalar data
Besides gridded data, it is also possible to write scalar data to a NetCDF file. Below is an
example that writes scalar data to the file "output\_scalar\_moselle.nc". For each NetCDF
variable a `name` (external variable name) and `parameter` (internal model parameter) is
required. A `reducer` can be specified to apply to the model output, see for more
information the following section [Output CSV section](@ref). When a `map` is provided to
extract data for certain locations (e.g. `gauges`) or areas (e.g. `subcatchment`), the
NetCDF location names are extracted from these maps. For a specific location (grid cell) a
`location` is required. Model parameters with the extra dimension `layer` for layered model
parameters of the vertical `sbm` concept require the specification of the layer (see also
example below). If multiple layers are desired, this can be specified in separate
`[[netcdf.variable]]` entries. In the section [Output CSV section](@ref), similar
functionality is available for CSV. For integration with Delft-FEWS, see also [Run from
Delft-FEWS](@ref), it is recommended to write scalar data to NetCDF format since the General
Adapter of Delft-FEWS can ingest this data format directly.

```toml
[netcdf]
path = "output_scalar_moselle.nc"

[[netcdf.variable]]
name = "Q"
map = "gauges"
parameter = "lateral.river.q"

[[netcdf.variable]]
coordinate.x = 6.255
coordinate.y = 50.012
name = "vwc_layer2_bycoord"
location = "vwc_bycoord"
parameter = "vertical.vwc"
layer = 2

[[netcdf.variable]]
location = "temp_byindex"
name = "temp_index"
index.x = 100
index.y = 264
parameter = "vertical.temperature"
```

## Output CSV section
Model output can also be written to CSV output. Below is an example that writes model output
to the file "output_moselle.csv". For each CSV column a `header` and `parameter` (internal
model parameter) is required. A `reducer` can be specified to apply to the model output,
with the following available reducers:

+ maximum
+ minimum
+ mean
+ median
+ first
+ last
+ only

with `only` as the default. To extract data for a specific location (grid cell), the `index`
of the vector, the coordinates `coordinate.x` and `coordinate.y`, or the x and y indices of
the 2D array (`index.x` and `index.y`) can be provided. Finally a `map` can be provided to
extract data for certain locations (e.g. `gauges`) or areas (e.g. `subcatchment`). In this
case a single entry can lead to multiple columns in the CSV file, which will be of the form
`header_id`, e.g. `Q_20`, for a gauge with integer ID 20. Model parameters with the extra
dimension `layer` for layered model parameters of the vertical `sbm` concept require the
specification of the layer (see also example below). If multiple layers are desired, this
can be specified in separate `[[csv.column]]` entries.

The double brackets in `[[csv.column]]` is TOML syntax to indicate that it is part of a
list. You may specify as many entries as you wish.

```toml
[csv]
path = "output_moselle.csv"

[[csv.column]]
header = "Q"
parameter = "lateral.river.q"
reducer = "maximum"

[[csv.column]]
header = "volume"
index = 1
parameter = "lateral.river.reservoir.volume"

[[csv.column]]
coordinate.x = 6.255
coordinate.y = 50.012
header = "temp_bycoord"
parameter = "vertical.temperature"

[[csv.column]]
coordinate.x = 6.255
coordinate.y = 50.012
header = "vwc_layer2_bycoord"
parameter = "vertical.vwc"
layer = 2

[[csv.column]]
header = "temp_byindex"
index.x = 100
index.y = 264
parameter = "vertical.temperature"

[[csv.column]]
header = "Q"
map = "gauges"
parameter = "lateral.river.q"

[[csv.column]]
header = "recharge"
map = "subcatchment"
parameter = "vertical.recharge"
reducer = "mean"
```

## Modify parameters
It is possible to modify model parameters and forcing through the TOML file. Two options to
modify input parameters are available:

- Set an input parameter (static) to an uniform value.
- Modify an input parameter (cyclic and static) or forcing variable through the use of a
  `scale` factor and `offset`.

To set for example the input parameter `cfmax` to an uniform value of 2.5:

```toml
[input.vertical]
water_holding_capacity = "WHC"
waterfrac = "WaterFrac"
cfmax.value = 2.5
```

To change for example the forcing variable `precipitation` with a `scale` factor of 1.5 and
an `offset` of 0.5:

```toml
[input.vertical.precipitation]
netcdf.variable.name = "P"
scale = 1.5
offset = 0.5
```

## Fixed forcing values
It is possible to set fixed values for forcing parameters through the TOML file. To set for
example `temperature` to a fixed value of 10 ``\degree``C, the complete `forcing` list is
required:

```toml
forcing = [
  "vertical.precipitation",
  "vertical.temperature",
  "vertical.potential_evaporation",
]

[input.vertical.temperature]
value = 10
```

Note that the mapping to the external netCDF variable listed under the section
`[input.vertical]` needs to be removed or commented out:

```toml
[input.vertical]
potential_evaporation = "PET" # forcing
# temperature = "TEMP" # forcing
precipitation = "P" # forcing
```