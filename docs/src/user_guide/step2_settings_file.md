# [Step 2: Preparing the settings file](@id config_toml)
A settings file is essential for wflow, as it contains information about the model
configuration, simulation period, input files, and parameters. The settings are provided in
a TOML file. The settings file is structured in several sections, which are explained below.
The file paths provided in this file are relative to the location of the TOML file, or to
`dir_input` and `dir_output` if they are specified.

## General time info
Time information is optional. When omitted, wflow will perform computations for each
timestamp in the forcing netCDF file, except for the first forcing timestamp, which is
considered equal to the initial conditions of the wflow model (state time). If you wish to
calculate a subset of this time range, or a different timestep, you can specify a
`starttime`, `endtime` and `timestepsecs`. The `starttime` is defined as the model state
time. In the TOML file settings below, the `starttime` is 2000-01-01T00:00:00 (state time)
and the first update (and output) of the wflow model is at 2000-01-02T00:00:00. The
`time_units` optional information is used by the `writer` of the model, for model output in
netCDF format. The `calendar` option allows you to calculate in one of the different [CF
conventions calendars](http://cfconventions.org/cf-conventions/cf-conventions.html#calendar)
provided by the [CFTime.jl package](https://juliageo.org/CFTime.jl/latest/), such as
`"360_day"`. This is useful if you want to calculate climate scenarios which are sometimes
provided in these alternative calendars.

```toml
calendar = "standard"                           # optional, this is the default value
starttime = 2000-01-01T00:00:00                 # optional, default from forcing netCDF
endtime = 2000-02-01T00:00:00                   # optional, default from forcing netCDF
time_units = "days since 1900-01-01 00:00:00"   # optional, this is the default value
timestepsecs = 86400                            # optional, default from forcing netCDF
dir_input = "data/input"                        # optional, default is the path of the TOML
dir_output = "data/output"                      # optional, default is the path of the TOML
```

## [Logging](@id logging_toml)
Wflow prints logging messages at various levels such as debug, info, and error. These
messages are sent to both the terminal and a log file. Note that logging to a file is only
part of the `Wflow.run(tomlpath::AbstractString)` method. If you want to debug an issue, it
can be helpful to set `loglevel = "debug"` in the TOML. To avoid flooding the screen, debug
messages are only sent to the log file. The following settings will affect the logging:

```toml
silent = false          # optional, default is "false"
loglevel = "debug"      # optional, default is "info"
path_log = "log.txt"    # optional, default is "log.txt"
fews_run = false        # optional, default value is false
```

`silent` avoids logging to the terminal, and only writes to the log file. `loglevel`
controls which levels are filtered out; for instance, the default setting `"info"` does not
print any debug-level messages. Note that for finer control, you can also pass an integer
log level. For details, see Julia's
[Logging](https://docs.julialang.org/en/v1/stdlib/Logging/#Log-event-structure)
documentation. `path_log` sets the desired output path for the log file. For information on
`fews_run`, see [Run from Delft-FEWS](@ref run_fews).

## Model section
Model-specific settings can be included in the model section of the TOML file.

```toml
[model]
type = "sbm"                        # one of ("sbm" or "sbm_gwf)
masswasting = false                 # include lateral snow transport in the model, default is false
snow = false                        # include snow modelling, default is false
reinit = true                       # cold (reinit = true) or warm state (reinit = false), default is true
reservoirs = false                  # include reservoir modelling, default is false
kin_wave_iteration = false          # enable kinematic wave iterations in the model, default is false
thicknesslayers = [100, 300, 800]   # specific SBM setting: for each soil layer, a thickness [mm] is specified
min_streamorder_river = 5           # minimum stream order to delineate subbasins for river domain, default is 6 (for multi-threading computing purposes)
min_streamorder_land = 4            # minimum stream order to delineate subbasins for land domain, default is 5 (for multi-threading computing purposes)

```

## State options
The `state` section in the TOML file provides information about the location of input and
output states of the model. This section is mostly relevant if the model needs to start with
a "warm" state (i.e. based on the results of a previous simulation). The example below shows
how to save the output states of the current simulation, so it can be used to initialize
another model in the future. Details on the settings required to start a model with a warm
state can be found in the [additional model options](@ref reinit). If it is not required to
store the outstates of the current simulation, the entire `state` section can be removed.

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
parameters files (in netCDF format). Forcing is applied to the vertical component of the
model, and needs to be mapped to the external netCDF variable name. `forcing` lists the
internal model forcing parameters, and these are mapped to the external netCDF variables
listed under the section `[input.vertical]`. It is possible to provide cyclic parameters to
the model (minimum time step of 1 day). In the example below, the internal
`vertical.leaf_area_index` model parameter is mapped to the external netCDF variable "LAI"
variable. Cyclic time inputs of parameters can be different (e.g., daily or monthly). The
`time` dimension name of these cylic input parameters in the model parameter netCDF file
should start with "time". If a model parameter is not mapped, a default value will be used,
if available.

```toml
[input]
# use "forcing-year-*.nc" if forcing files are split in time
path_forcing = "forcing-moselle.nc"    # Location of the forcing data
path_static = "staticmaps-moselle.nc"  # Location of the static data

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

[input.vertical]    # Map internal model variable/parameter names to variable names in the netCDF files
altitude = "wflow_dem"
c = "c"
cf_soil = "cf_soil"
cfmax = "Cfmax"
e_r = "EoverR"
infiltcappath = "InfiltCapPath"
infiltcapsoil = "InfiltCapSoil"
kext = "Kext"
kv_0 = "KsatVer"
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
theta_r = "thetaR"
theta_s = "thetaS"

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

## Output netCDF section

### Grid data
This optional section of the TOML file specifies the output netCDF file for writing gridded
model output. It includes a mapping between internal model parameter components and external
netCDF variables.

To limit the size of the resulting netCDF file, file compression can be enabled. Compression
increases computational time but can significantly reduce the size of the netCDF file. Set
the `compressionlevel` variable to a value between `0` and `9`. A setting of `0` means no
compression, while values between 1 and 9 indicate increasing levels of compression (1:
least compression, minimal run-time impact, 9: highest compression, maximum run-time
impact). If file size is a concern, we recommend using a value of `1`, as higher compression
levels generally have a limited effect on file size.

```toml
[output]
path = "output_moselle.nc"         # Location of the output file
compressionlevel = 1               # Compression level (default 0)

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
In addition to gridded data, scalar data can also be written to a netCDF file. Below is an
example that shows how to write scalar data to the file "output\_scalar\_moselle.nc". For
each netCDF variable, a `name` (external variable name) and a `parameter` (internal model
parameter) are required. A `reducer` can be specified to apply to the model output. See more
details in the [Output CSV section](@ref) section. If a `map` is provided to extract data
for specific locations (e.g. `gauges`) or areas (e.g. `subcatchment`), the netCDF location
names are extracted from these maps. For a specific location (grid cell) a `location` is
required. For layered model parameters and variables that have an extra `layer` dimension
and are part of the vertical `sbm` concept, an internal layer index can be specified (an
example is provided below). If multiple layers are desired, this can be specified in
separate `[[netcdf.variable]]` entries. Note that the additional dimension should be
specified when wflow is integrated with Delft-FEWS, for netCDF scalar data an extra
dimension is not allowed by the `importNetcdfActivity` of the Delft-FEWS General Adapter. In
the section [Output CSV section](@ref), similar functionality is available for CSV. For
integration with Delft-FEWS, it is recommended to write scalar data to netCDF format, as the
General Adapter of Delft-FEWS can directly ingest data from netCDF files. For more
information, see [Run from Delft-FEWS](@ref run_fews). 

```toml
[netcdf]
path = "output_scalar_moselle.nc"  # Location of the results

[[netcdf.variable]] # Extract the values of lateral.river.q using the gauges map, and assign it with the name 'Q' as a variable in the netCDF file
name = "Q"
map = "gauges"
parameter = "lateral.river.q"

[[netcdf.variable]] # Using coordinates to extract temperature
coordinate.x = 6.255
coordinate.y = 50.012
name = "vwc_layer2_bycoord"
location = "vwc_bycoord"
parameter = "vertical.vwc"
layer = 2

[[netcdf.variable]] # Using indices to extract temperature
location = "temp_byindex"
name = "temp_index"
index.x = 100
index.y = 264
parameter = "vertical.temperature"
```

## Output CSV section
Model output can also be written to a CSV file. Below is an example that writes model output
to the file "output_moselle.csv". For each CSV column, a `header` and `parameter` (internal
model parameter) are required. A `reducer` can be specified to apply to the model output,
with the following available reducers:

+ maximum
+ minimum
+ mean
+ median
+ sum
+ first
+ last
+ only

with `only` as the default. To extract data for a specific location (grid cell), the `index`
of the vector, the coordinates `coordinate.x` and `coordinate.y`, or the x and y indices of
the 2D array (`index.x` and `index.y`) can be provided. Additionally, a `map` can be
provided to extract data for certain locations (e.g. `gauges`) or areas (e.g.
`subcatchment`). In this case, a single entry can lead to multiple columns in the CSV file,
which will be of the form `header_id`, e.g. `Q_20`, for a gauge with integer ID 20. For
layered model parameters and variables that have an extra dimension `layer` and are part of
the vertical `sbm` concept an internal layer index (see also example below) should be
specified. If multiple layers or classes are desired, this can be specified in separate
`[[csv.column]]` entries.

The double brackets in `[[csv.column]]` follow TOML syntax, indicating that it is part of a
list. You can specify as many entries as you want.

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

- Set an input parameter (static) to a uniform value.
- Modify an input parameter (cyclic and static) or forcing variable using a `scale` factor
  and `offset`.

For example, to set the input parameter `cfmax` to an uniform value of 2.5:

```toml
[input.vertical]
water_holding_capacity = "WHC"
waterfrac = "WaterFrac"
cfmax.value = 2.5
```

For input parameters with an extra dimension (e.g. `layer` or `classes`), one uniform value
can be provided or a list of values that matches the length of the additional dimension. For
example, a list of values can be provided for input parameter `c` as follows:

```toml
[input.vertical]
water_holding_capacity = "WHC"
waterfrac = "WaterFrac"
c.value = [10.5, 11.25, 9.5, 7.0]
```

To change the forcing variable `precipitation` with a `scale` factor of 1.5 and
an `offset` of 0.5:

```toml
[input.vertical.precipitation]
netcdf.variable.name = "P"
scale = 1.5
offset = 0.5
```

For input parameters with an extra dimension, it is also possible to modify multiple indices
simultaneously with different `scale` and `offset` values. In the example below, the
external netCDF variable `c` is modified at `layer` index 1 and 2, with a `scale` factor of
2.0 and 1.5 respectively, and an `offset` of 0.0 for both indices:

```toml
[input.vertical.c]
netcdf.variable.name = "c"
scale = [2.0, 1.5]
offset = [0.0, 0.0]
layer = [1, 2]
```

## Fixed forcing values
It is possible to set fixed values for forcing parameters through the TOML file. For
example, to set `temperature` to a fixed value of 10 ``\degree``C, the complete `forcing`
list is required:

```toml
forcing = [
  "vertical.precipitation",
  "vertical.temperature",
  "vertical.potential_evaporation",
]

[input.vertical.temperature]
value = 10
```

Note that the mapping to the external netCDF variable listed under the `[input.vertical]`
section needs to be removed or commented out:

```toml
[input.vertical]
potential_evaporation = "PET" # forcing
# temperature = "TEMP" # forcing
precipitation = "P" # forcing
```
