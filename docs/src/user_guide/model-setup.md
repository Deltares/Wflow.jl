# Building a model from scratch

## Data requirements
The actual data requirements depend on the application of the Model and the Model type. Both
forcing and static data should be provided in netCDF format, with the same grid definition
for forcing and static data. The only exception is storage and rating curves for lakes, that
should be provided in CSV format, see also [Additional settings](@ref).

* Forcing data:
  - Precipitation
  - Potential evapotranspiration
  - Temperature (optional, only needed for snow and glacier modelling)

The requirements for static data (including model parameters) depend on the Model type. The
following data is required for all Model types, but not directly part of a Model component:

+ flow direction data (D8)
+ river map (location of the river)
+ sub-catchment map (model domain)

For the flow direction (D8) data, the PCRaster `ldd` convention is used, see also [PCRaster
ldd](https://pcraster.geo.uu.nl/pcraster/4.3.1/documentation/pcraster_manual/sphinx/secdatbase.html#ldd-data-type).
An approach to generate `ldd` data is to make use of the Python package
[pyflwdir](https://github.com/Deltares/pyflwdir):

+ to [upscale existing flow direction
  data](https://deltares.github.io/pyflwdir/latest/_examples/upscaling.html) as the 3 arcsec MERIT
  Hydro data (Yamazaki et al., 2019)
+ or to [derive flow directions from elevation
  data](https://deltares.github.io/pyflwdir/latest/_examples/from_dem.html),

see also Eilander et al. (2021) for more information.
Pyflwdir is also used by the [hydroMT](@ref) Python package described in the next paragraph.
Another approach to generate `ldd` data is to make use of PCRaster functionality, see for
example
[lddcreate](https://pcraster.geo.uu.nl/pcraster/4.3.1/documentation/pcraster_manual/sphinx/op_lddcreate.html).

Optionally, but also not directly part of a model component are `gauge` locations, that are
used to extract gridded data from certain locations.

The different supported model configurations are described in the section [Model
configurations](@ref). Wflow\_sbm models have the vertical concept [SBM](@ref vert_sbm) in
common and input parameters for this component are described in the [SBM](@ref params_sbm)
section of Model parameters. For wflow\_sbm models there are two ways to include subsurface
flow:

1. The kinematic wave approach (see section [Subsurface flow routing](@ref)) as part of the
   `sbm` model type. Parameters that are part of this component are described in the
   [Lateral subsurface flow](@ref params_ssf) section of Model parameters. Input parameters
   for this component are derived from the SBM vertical concept and the land slope. One
   external parameter [`ksathorfrac`](@ref params_ssf) is used to calculate the horizontal
   hydraulic conductivity at the soil surface `kh_0`.
2. Groundwater flow (see section [Groundwater flow component](@ref lateral_gwf)) as part of
   the `sbm_gwf` model type. For the unconfined aquifer the input parameters are described
   in the section [Unconfined aquifer](@ref) of Model parameters. The bottom (`bottom`) of
   the groundwater layer is derived from from the `soilthickness` [mm] parameter of `SBM`
   and the provided surface elevation `altitude` [m] as part of the static input. The `area`
   parameter is derived from the model grid. Parameters that are part of the boundary
   conditions of the unconfined aquifer are listed under [Constant Head](@ref) and [Boundary
   conditions](@ref) of the Model parameters section.

Most hydrological model configurations make use of the kinematic wave surface routing (river
flow, overland flow or both) and input data required for the river and overland flow
components is described in [Surface flow](@ref).  There is also the option to use the local
inertial model as part of the wflow\_sbm models (model types `sbm` and `sbm_gwf`):
+ for river flow, see also the [Local inertial river and floodplain](@ref
  config_sbm_gwf_lie_river) model.
+ for 1D river flow and 2D overland flow combined, see also the [Local inertial river (1D)
  and land (2D)](@ref config_sbm_gwf_lie_river_land) model.

Input parameters for this approach are described in [River flow (local inertial)](@ref
local-inertial_river_params), including the optional 1D [floodplain schematization](@ref
local-inertial_floodplain_params), and [Overland flow (local inertial)](@ref
local-inertial_land_params) of the Model parameters section.

Reservoirs or lakes can be part of the kinematic wave or local inertial model for river flow
and input parameters are described in [Reservoirs](@ref reservoir_params) and [Lakes](@ref
lake_params).

The [wflow\_hbv](@ref config_hbv) model configuration consists besides the river and
overland flow components of the [HBV](@ref vert_hbv) vertical concept. Input parameters for
this component are described in the [HBV](@ref params_hbv) section of Model parameters. For
the river and overland flow components the kinematic wave approach is used.

The [wflow\_flextopo](@ref config_flextopo) model configuration consists besides the river
and overland flow components of the [FLEXTopo](@ref vert_flextopo) vertical concept. Input
parameters for this component are described in the [FLEXTopo](@ref params_flextopo) section
of Model parameters. For the river and overland flow components the kinematic wave approach
is used.

The [wflow\_sediment](@ref config_sediment) model configuration consists of the vertical
[Soil Erosion](@ref) concept and the input parameters for this concept are described in the
[Sediment](@ref params_sediment) section of the Model parameters. The parameters of the
lateral [Sediment Flux in overland flow](@ref) concept are described in the [Overland
flow](@ref) section of the Model parameters. Parameters of this component are not directly
set by data from static input. The input parameters of the lateral concept [River Sediment
Model](@ref) are listed in [River flow](@ref) of the Model parameters section.

The Model parameters section lists all the parameters per Model component and these Tables
can also be used to check which parameters can be part of the output, see also [Output
netCDF section](@ref) and [Output CSV section](@ref).

Example models can be found in the [Example model section](@ref sample_data).

## hydroMT
[hydroMT](https://github.com/Deltares/hydromt) is a Python package, developed by Deltares,
to build and analyze hydro models. It provides a generic model api with attributes to
access the model schematization, (dynamic) forcing data, results and states.

For the following wflow\_sbm model (modeltype `sbm`) configurations:
  - [wflow\_sbm + kinematic wave routing](@ref config_sbm)
  - [wflow\_sbm + local inertial river and floodplain](@ref config_sbm_gwf_lie_river)
  - [wflow\_sbm + local inertial river (1D) and land (2D)](@ref config_sbm_gwf_lie_river_land)
and the [wflow\_sediment](@ref config_sediment) model configuration, the wflow plugin
[hydroMT-wflow](https://github.com/Deltares/hydromt_wflow) of hydroMT can be used to build
and analyze these wflow models in an automated way.

To learn more about the wflow plugin of this Python package, we refer to the [hydroMT-wflow
documentation](https://deltares.github.io/hydromt_wflow/latest/index.html).

To inspect or modify (for example in QGIS) the netCDF static data of these wflow models it
is convenient to export the maps to a raster format. This can be done as part of the
hydroMT-wflow plugin, see also the following [example]
(https://deltares.github.io/hydromt_wflow/latest/examples/examples/convert_staticmaps_to_mapstack.html).
It is also possible to create again the netCDF static data file based on the modified raster
map stack.


## References
+ Yamazaki, D., Ikeshima, D., Sosa, J., Bates, P. D., Allen, G. H. and Pavelsky, T. M.:
  MERIT Hydro: A high‐resolution global hydrography map based on latest topography datasets,
  Water Resour. Res., 2019WR024873, doi:10.1029/2019WR024873, 2019.
+ Eilander, D., van Verseveld, W., Yamazaki, D., Weerts, A., Winsemius, H. C., and Ward, P.
  J.: A hydrography upscaling method for scale-invariant parametrization of distributed
  hydrological models, Hydrol. Earth Syst. Sci., 25, 5287–5313,
  <https://doi.org/10.5194/hess-25-5287-2021>, 2021.
