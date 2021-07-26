# Building a model

## Data requirements
The actual data requirements depend on the application of the model and the model type. Both
forcing and static data should be provided in netCDF format, with the same grid definition
for forcing and static data. The only exception is storage and rating curves for lakes, that
should be provided in CSV format, see also [Additional settings](@ref).

* Forcing data:
  - Precipitation
  - Potential evapotranspiration
  - Temperature (optional, only needed for snow and glacier modelling)

The requirements for static data (including model parameters) depend on the model type. The 
following data is required for all model types, but not directly part of a model component:

+ flow direction data (D8)
+ river map (location of the river)
+ sub-catchment map (model domain)

For the flow direction (D8) data, the PCRaster `ldd` convention is used, see also [PCRaster
ldd](https://pcraster.geo.uu.nl/pcraster/4.3.1/documentation/pcraster_manual/sphinx/secdatbase.html#ldd-data-type).
An approach to generate `ldd` data is to make use of the Python package
[pyflwdir](https://deltares.gitlab.io/wflow/pyflwdir/) to [upscale existing flow direction
data](https://deltares.gitlab.io/wflow/pyflwdir/flwdir.html#Flow-direction-upscaling) as the
3 arcsec MERIT Hydro data (Yamazaki et al., 2019), see also Eilander et al. (2020) for more
information. Pyflwdir is also used by the [hydroMT](@ref) Python package described in the
next paragraph. Another approach to generate `ldd` data is to make use of PCRaster
functionality, see for example
[lddcreate](https://pcraster.geo.uu.nl/pcraster/4.3.1/documentation/pcraster_manual/sphinx/op_lddcreate.html).

Optionally, but also not directly part of a model component are `gauge` locations, that are
used to extract gridded data from certain locations.

The model types that make use of the kinematic wave:
+ [SBM + Kinematic wave](@ref)
+ [SBM + Groundwater flow](@ref)
+ [HBV model](@ref)

require for the river and overland flow components input data that is described in [Surface
flow](@ref). Reservoirs or lakes can be part of the kinematic wave (optional) and input
parameters are described in [Reservoirs](@ref reservoir_params) and [Lakes](@ref
lake_params).

Besides the river and overland flow components the [SBM + Kinematic wave](@ref) model
consists of the [SBM vertical concept](@ref) and input parameters for this component are
described in the [SBM](@ref) section of Model parameters. Finally, the [SBM + Kinematic
wave](@ref) model includes the lateral component [Subsurface flow routing](@ref) and
parameters that are part of this component are described in the [Lateral subsurface
flow](@ref) section of Model parameters.

## hydroMT
[hydroMT](https://github.com/Deltares/hydromt) is a Python package, developed by Deltares,
to build and analysis hydro models. It provides a generic model api with attributes to
access the model schematization, (dynamic) forcing data, results and states.

For the following Wflow models:
  - [SBM + Kinematic wave](@ref) 
  - [wflow_sediment](@ref)

the Wflow plugin [hydroMT-wflow](https://github.com/Deltares/hydromt_wflow) of hydroMT can
be used to build and analyse these Wflow model types in an automated way.

To learn more about the Wflow plugin of this Python package, we refer to the [hydroMT-wflow
documentation](https://deltares.github.io/hydromt_wflow/latest/index.html).

To inspect or modify (for example in QGIS) the netCDF static data of these Wflow models it
is convenient to export the maps to a raster format. This can be done as part of the
hydroMT-wflow plugin, see also the following [example]
(https://deltares.github.io/hydromt_wflow/latest/examples/examples/convert_staticmaps_to_gtiff.html).
It is also possible to create again the netCDF static data file based on the modified raster
map stack.


## References
+ Yamazaki, D., Ikeshima, D., Sosa, J., Bates, P. D., Allen, G. H. and Pavelsky, T. M.:
  MERIT Hydro: A high‐resolution global hydrography map based on latest topography datasets,
  Water Resour. Res., 2019WR024873, doi:10.1029/2019WR024873, 2019.
+ Eilander, D., van Verseveld, W., Yamazaki, D., Weerts, A., Winsemius, H. C., and Ward, P.
  J.: A hydrography upscaling method for scale invariant parametrization of distributed
  hydrological models, Hydrol. Earth Syst. Sci. Discuss. [preprint],
  <https://doi.org/10.5194/hess-2020-582>, in review, 2020. 