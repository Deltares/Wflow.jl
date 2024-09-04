# Additional wflow options

## [Starting the model with "warm" states](@id reinit)

The `state` section in the TOML file provides information on the input file if the model is
initialized with a warm state (`path_input`) and to what file the states are written at the
end of the model run (`path_output`). Please note that the model setting `reinit` needs to
be set to `false` in order to initialize the model with states from the file located at
`path_input`. A mapping between external state names and internal model states is required.
This information is specified for each model component, the `vertical` model and `lateral`
model components. In the example below the `vertical` component represents the SBM concept,
and for the `lateral` components there is a `river` (including optional `reservoir`, `lake`
and `floodplain` components), `land` and `subsurface` domain. The internal model states are
listed on the left side, and the external state names are listed on the right side. Note
that `path_input` is only required when `reinit` is set to false. `path_output` is optional,
an output state file is only written when it is defined. If neither is set, the entire
`state` section can be left out.

```toml
[model]
reinit = false # cold (reinit = true) or warm state (reinit = false), default is true

[state]
path_input = "data/instates-moselle.nc"     # Location of the file with the input states
path_output = "data/outstates-moselle.nc"   # Output location of the states after the model run

[state.vertical]
satwaterdepth = "satwaterdepth"
snow = "snow"
tsoil = "tsoil"
ustorelayerdepth = "ustorelayerdepth"
canopystorage = "canopystorage"
snowwater = "snowwater"
glacierstore ="glacierstore"

[state.lateral.river]
q = "q_river"
h = "h_river"
h_av = "h_av_river"

[state.lateral.river.floodplain] 
q = "q_floodplain"
h = "h_floodplain"

[state.lateral.river.reservoir]
volume = "volume_reservoir"

[state.lateral.river.lake]
waterlevel = "waterlevel_lake"

[state.lateral.subsurface]
ssf = "ssf"

[state.lateral.land]
q = "q_land"
h = "h_land"
h_av = "h_av_land"
```

## Enabling snow and glacier processes

```toml
[model]
snow = true
masswasting = true
glacier = true

[input.vertical]
tt = "TT"
tti = "TTI"
ttm = "TTM"
water_holding_capacity = "WHC"
glacierstore = "wflow_glacierstore"
glacierfrac = "wflow_glacierfrac"
g_cfmax = "G_Cfmax"
g_tt = "G_TT"
g_sifrac = "G_SIfrac"
```

## Enabling reservoirs

```toml
[model]
reservoirs = true

[input.lateral.river.reservoir]
area = "ResSimpleArea"
areas = "wflow_reservoirareas"
demand = "ResDemand"
locs = "wflow_reservoirlocs"
maxrelease = "ResMaxRelease"
maxvolume = "ResMaxVolume"
targetfullfrac = "ResTargetFullFrac"
targetminfrac = "ResTargetMinFrac"

[state.lateral.river.reservoir]
volume = "volume_reservoir"
```

## Enabling lakes

```toml
[model]
lakes = true

[input.lateral.river.lake]
area = "lake_area"
areas = "wflow_lakeareas"
b = "lake_b"
e = "lake_e"
locs = "wflow_lakelocs"
outflowfunc = "lake_outflowfunc"
storfunc  = "lake_storfunc"
threshold  = "lake_threshold"
waterlevel = "lake_waterlevel"

[state.lateral.river.lake]
waterlevel = "waterlevel_lake"
```

## Enabling Floodplain routing
As part of the local inertial model for river flow.

```toml
[model]
floodplain_1d = true

[input.lateral.river.floodplain]
volume = "floodplain_volume"
n = "floodplain_n"

[state.lateral.river.floodplain] 
q = "q_floodplain"
h = "h_floodplain"
```

## Enabling water demand and allocation
The model types `sbm` and `sbm_gwf` support water demand and allocation computations, in
combination with the kinematic wave and local inertial runoff routing scheme for river and
overland flow.

```toml
# example of water demand and allocation input parameters as cyclic data
[input]
cyclic = ["vertical.domestic.demand_gross", "vertical.domestic.demand_net", 
"vertical.industry.demand_gross", "vertical.industry.demand_net", 
"vertical.livestock.demand_gross", "vertical.livestock.demand_net", 
"vertical.paddy.irrigation_trigger", "vertical.nonpaddy.irrigation_trigger",]

[model.water_demand]
domestic = true     # optional, default is "false"
industry = true     # optional, default is "false"
livestock = true    # optional, default is "false"
paddy = true        # optional, default is "false"
nonpaddy = true     # optional, default is "false"

[input.vertical.allocation]
areas = "allocation_areas"
frac_sw_used = "SurfaceWaterFrac"

[input.vertical.domestic]
demand_gross = "dom_gross"
demand_net = "dom_net"

[input.vertical.industry]
demand_gross = "ind_gross"
demand_net = "ind_net"

[input.vertical.livestock]
demand_gross = "lsk_gross"
demand_net = "lsk_net"

[input.vertical.paddy]
irrigation_areas = "paddy_irrigation_areas"
irrigation_trigger = "irrigation_trigger"

[input.vertical.nonpaddy]
irrigation_areas = "nonpaddy_irrigation_areas"
irrigation_trigger = "irrigation_trigger"

# required if paddy is set to "true"
[state.vertical.paddy]
h = "h_paddy"
```

## [Using multithreading] (@id multi_threading)

### Using wflow in Julia

Wflow supports multi-threading execution of the wflow\_sbm model that uses the kinematic
wave approach for river, overland and lateral subsurface flow. Both the vertical SBM concept
and the kinematic wave components of this model can run on multiple threads. The optional
[local inertial model for river flow](@ref config_sbm_gwf_lie_river_land) and the optional
[local inertial model for river (1D) and land (2D)](@ref config_sbm_gwf_lie_river_land),
both part of wflow\_sbm, can also run on multiple threads. The threading functionality for
the kinematic wave may also be useful for the wflow\_sbm model [SBM + Groundwater flow](@ref
config_sbm_gwf). The multi-threading functionality in wflow is considered experimental, see
also the following [issue](https://github.com/Deltares/Wflow.jl/issues/139), where an error
was not thrown running code multi-threaded. Because of this we advise to start with running
a wflow model single-threaded (for example during the testing phase of setting up an new
wflow model).

For information on how to start Julia with multiple threads we refer to [How to start Julia
with multiple
threads](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads).

Additionally, when running Julia + wflow via the command line (note that this is different
from the `wflow_cli`), it is possible to define the number of threads via the `-t` flag.
An example where we start Julia with three threads:

```
julia -t 3 -e 'using Wflow; Wflow.run()' path/to/config.toml
```

### [Using the command line interface](@id cli_multi_threading)

As explained above, we need to start julia with multiple threads to make use of this
speedup. For `wflow_cli`, the only way to do this is by setting the `JULIA_NUM_THREADS`
environment variable, as explained in [these julia
docs](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads).

When a model run starts, among the run information the number of threads that are used is
printed, so `nthreads() = 4` means 4 threads are used, because `JULIA_NUM_THREADS` has been
set to 4.

## Using the Basic Model Interface

### Introduction
The [Community Surface Dynamics Modeling System](https://csdms.colorado.edu/wiki/Main_Page)
(CSMDS) has developed the Basic Model Interface (BMI). BMI consists of a set of standard
control and query functions that can be added by a developer to the model code and makes a
model both easier to learn and easier to couple with other software elements.

For more information see also: <http://csdms.colorado.edu/wiki/BMI_Description>

CSDMS provides specifications for the languages C, C++, Fortran and Python. Wflow, written
in the [Julia programming language](https://julialang.org/), makes use of the following
[Julia specification](https://github.com/Deltares/BasicModelInterface.jl), based on BMI 2.0
version.

For the BMI implementation of wflow all grids are defined as [unstructured
grids](https://bmi-spec.readthedocs.io/en/latest/model_grids.html#unstructured-grids),
including the special cases `scalar` and `points`. While the input (forcing and model
parameters) is structured (uniform rectilinear), internally wflow works with one dimensional
arrays based on the active grid cells of the 2D model domain.

### Configuration
The variables that wflow can exchange through BMI are based on the different model
components and these components should be listed under the `API` section of the TOML
configuration file of the model type. Below an example of this `API` section, that lists the
`vertical` component and different `lateral` components:

```toml
[API]
components = [
  "vertical",
  "lateral.subsurface",
  "lateral.land",
  "lateral.river",
  "lateral.river.reservoir"
]
```

See also:
```@docs
Wflow.BMI.initialize
Wflow.BMI.get_input_var_names
```

Variables with a third dimension, for example `layer` as part of the vertical `SBM` concept,
are exposed as two-dimensional grids through the wflow BMI implementation. For these
variables the index of this third dimension is required, by adding `[k]` to the variable
name (`k` refers to the index of the third dimension). For example, the variable
`vertical.vwc[1]` refers to the first soil layer of the vertical `SBM` concept.

### Couple to a groundwater model
For the coupling of wflow\_sbm (SBM + kinematic wave) with a groundwater model (e.g.
MODFLOW) it is possible to run:
- wflow\_sbm in parts from the BMI, and
- to switch off the lateral subsurface flow component of wflow\_sbm.

The lateral subsurface component of wflow\_sbm is not initialized by wflow when the
`[input.lateral.subsurface]` part of the TOML file is not included. Then from the BMI it is
possible to run first the recharge part of SBM:

```julia
model = BMI.update(model, run="sbm_until_recharge")
```
and to exchange recharge and for example river waterlevels to the groundwater model. After
the groundwater model update, and the exchange of groundwater head and for example drain and
river flux to wflow\_sbm, the SBM part that mainly determines exfiltration of water from the
unsaturated store, and the kinematic wave for river - and overland flow can be run as
follows:

```julia
model = BMI.update(model, run="sbm_after_subsurfaceflow")
```

See also:
```@docs
Wflow.BMI.update
```

## [Run from Delft-FEWS](@id run_fews)

Wflow integrates easily as part of an operational system by linking to the
[Delft-FEWS](https://oss.deltares.nl/web/delft-fews/) platform. Delft-FEWS integrates data
and models, and is for example used in many active flood forecasting systems around the
world.

This can be done without a model adapter that provides the interface between Delft-FEWS and
an external model (or module). This is possible because time information in the TOML
configuration file is optional and Delft-FEWS can import and export netCDF files. When time
information is left out from the TOML configuration file, the `starttime`, `endtime` and
`timestepsecs` (timestep) of the run is extracted from the netCDF forcing file by wflow.

To indicate that a wflow model runs from Delft-FEWS, the following setting needs to be
specified in the main section of the TOML configuration file:

```toml
fews_run = true  # optional, default value is false
```

This ensures that wflow offsets the time handling, to meet the expectations of Delft-FEWS.

It also uses a different format for the log file such that each log message takes up only
one line. That meets the [General Adapter
logFile](https://publicwiki.deltares.nl/display/FEWSDOC/05+General+Adapter+Module#id-05GeneralAdapterModule-logFile)
expectations, which then can get parsed with these Delft-FEWS log parsing settings:

```
<logFile>
    <file>log.txt</file>
    <errorLinePattern >* [Error] *</errorLinePattern >
    <warningLinePattern>* [Warn] *</warningLinePattern>
    <infoLinePattern>* [Info] *</infoLinePattern>
    <debugLinePattern >* [Debug] *</debugLinePattern >
</logFile>
```

## Run wflow as a ZMQ Server
It is possible to run wflow as a ZMQ Server, for example for the coupling to the
[OpenDA](https://openda.org/) software for data-assimilation. The code for the wflow ZMQ
Server is not part of the Wflow.jl package, and is located
[here](https://github.com/Deltares/Wflow.jl/tree/master/server).
