# [Step 2: Preparing the settings file](@id config_toml)
A settings file is essential for wflow, as it contains information on the model configuration, simulation period, input files, and parameters. The settings are provided in a TOML file. The settings file is structured in several sections, which are explained below. The filepaths that are provided in this file are relative to the location of the TOML file. 

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
starttime = 2000-01-01T00:00:00                 # optional, default from forcing NetCDF
endtime = 2000-02-01T00:00:00                   # optional, default from forcing NetCDF
time_units = "days since 1900-01-01 00:00:00"   # optional, this is default value
timestepsecs = 86400                            # optional, default from forcing NetCDF
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
masswasting = false                 # include lateral snow transport in the model, default is false
snow = false                        # include snow modelling, default is false
reinit = true                       # cold (reinit = true) or warm state (reinit = false), default is true
reservoirs = false                  # include reservoir modelling, default is false
kin_wave_iteration = false          # enable kinematic wave iterations in the model, default is false
thicknesslayers = [100, 300, 800]   # specific SBM setting: for each soil layer a thickness [mm] is specified
min_streamorder = 4                 # minimum stream order to delineate subbasins, default is 4 (for multi-threading computing purposes)
```

## State options
The `state` section in the TOML file provides information about the location of input and output states of the model. This section is mostly relevant if the model needs to be started with a "warm" state (i.e. based on the results of a previous simulation). The example below shows how to save the output states of the current simulation, so it can be used to initialize another model in the future. Details on the settings required to start a model with a warm state can be found in the [additional model options](@ref reinit). If it is not required to store the outstates of the current simulation, the entire `state` section can be removed.

```toml
[state]
path_output = "data/outstates-moselle.nc"   # Output location of the states after the model run
```

## Input section
The `input` section of the TOML file contains information about the input forcing and model
parameters files (netCDF format). Forcing is applied to the vertical component of the model,
and needs to be mapped to the external netCDF variable name. `forcing` lists the internal
model forcing parameters, and these are mapped to the external netCDF variables listed under
the section `[input.vertical]`. It is possible to provide cyclic parameters to the model. In
the example below this is done for the internal `vertical.leaf_area_index` model parameter,
that is linked to the external netCDF variable "LAI" variable. 
If a model parameter is not mapped, a default value will be used if available.

```toml
[input]
path_forcing = "data/forcing-moselle.nc"    # Location of the forcing data
path_static = "data/staticmaps-moselle.nc"  # Location of the static data

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

[input.vertical]    # Map internal model variable/parameter names to names of the variables in the NetCDF files
altitude = "wflow_dem" 
c = "c" 
cf_soil = "cf_soil" 
cfmax = "Cfmax" 
e_r = "EoverR" 
infiltcappath = "InfiltCapPath" 
infiltcapsoil = "InfiltCapSoil" 
kext = "Kext" 
"kv₀" = "KsatVer" 
leaf_area_index = "LAI"             # Cyclic variable
m = "M" 
maxleakage = "MaxLeakage" 
pathfrac = "PathFrac" 
potential_evaporation = "PET"       # Forcing variable
precipitation = "P"                 # Forcing variable
rootdistpar = "rootdistpar" 
rootingdepth = "RootingDepth" 
soilminthickness = "SoilMinThickness" 
soilthickness = "SoilThickness" 
specific_leaf = "Sl" 
storage_wood = "Swood" 
temperature = "TEMP"                # Forcing variable
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
path = "data/output_moselle.nc"         # Location of the output file

[output.vertical]   # Mapping of names between internal model components and external netCDF variables
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
`location` is required. For layered model parameters and variables that have an extra
dimension `layer` and are part of the vertical `sbm` concept it is possible to specify an
internal layer index (see also example below). If multiple layers are desired, this can be
specified in separate `[[netcdf.variable]]` entries. Note that the specification of the
layer is not optional when Wflow is integrated with Delft-FEWS, for NetCDF scalar data an
extra dimension is not allowed by the `importNetcdfActivity` of the Delft-FEWS General
Adapter. In the section [Output CSV section](@ref), similar functionality is available for
CSV. For integration with Delft-FEWS, see also [Run from Delft-FEWS](@ref), it is
recommended to write scalar data to NetCDF format since the General Adapter of Delft-FEWS
can ingest this data format directly.

```toml
[netcdf]
path = "data/output_scalar_moselle.nc"  # Location of the results

[[netcdf.variable]] # Extract the values of lateral.river.q using the gauges map, and assigning it with the name 'Q' as variable to the netCDF
name = "Q"
map = "gauges"
parameter = "lateral.river.q"

[[netcdf.variable]] # Using coordinates to extract the temperature
coordinate.x = 6.255
coordinate.y = 50.012
name = "vwc_layer2_bycoord"
location = "vwc_bycoord"
parameter = "vertical.vwc"
layer = 2

[[netcdf.variable]] # Using indices to extract the temperature
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
`header_id`, e.g. `Q_20`, for a gauge with integer ID 20. For layered model parameters and
variables that have an extra dimension `layer` and are part of the vertical `sbm` concept an
internal layer index (see also example below) should be specified. If multiple layers are
desired, this can be specified in separate `[[csv.column]]` entries.

The double brackets in `[[csv.column]]` is TOML syntax to indicate that it is part of a
list. You may specify as many entries as you wish.

```toml
[csv]
path = "data/output_moselle.csv"

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

