# Building a model

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
[pyflwdir](https://deltares.gitlab.io/wflow/pyflwdir/) to [upscale existing flow direction
data](https://deltares.gitlab.io/wflow/pyflwdir/flwdir.html#Flow-direction-upscaling) as the
3 arcsec MERIT Hydro data (Yamazaki et al., 2019), see also Eilander et al. (2020) for more
information. Pyflwdir is also used by the [hydroMT](@ref) Python package described in the
next paragraph. Another approach to generate `ldd` data is to make use of PCRaster
functionality, see for example
[lddcreate](https://pcraster.geo.uu.nl/pcraster/4.3.1/documentation/pcraster_manual/sphinx/op_lddcreate.html).

Optionally, but also not directly part of a model component are `gauge` locations, that are
used to extract gridded data from certain locations.

The following Model types make use of the kinematic wave:
+ [SBM + Kinematic wave](@ref)
+ [SBM + Groundwater flow](@ref)
+ [HBV model](@ref)

and require for the river and overland flow components input data that is described in [Surface
flow](@ref). Reservoirs or lakes can be part of the kinematic wave (optional) and input
parameters are described in [Reservoirs](@ref reservoir_params) and [Lakes](@ref
lake_params).

Besides the river and overland flow components the SBM + Kinematic wave model consists of
the [SBM vertical concept](@ref) and input parameters for this component are described in
the [SBM](@ref) section of Model parameters. Finally, the SBM + Kinematic wave model
includes the lateral component [Subsurface flow routing](@ref) and parameters that are part
of this component are described in the [Lateral subsurface flow](@ref) section of Model
parameters. Input parameters for this component of the SBM + Kinematic wave model are
derived from the SBM vertical concept and the land slope. One external parameter `khfrac` is
used to calculate the horizontal hydraulic conductivity at the soil surface `kh₀`.

The HBV model consists besides the river and overland flow components of the [HBV vertical
concept](@ref). Input parameters for this component are described in the [HBV](@ref) section
of Model parameters.

The SBM + Groundwater flow includes besides the river and overland flow components and the
vertical SBM concept, the lateral [Groundwater flow component](@ref lateral_gwf). For the
unconfined aquifer the input parameters are described in the section [Unconfined
aquifer](@ref) of Model parameters. The bottom (`bottom`) of the groundwater layer is
derived from from the `soilthickness` [mm] parameter of `SBM` and the provided surface
elevation `altitude` [m] as part of the static input. The `area` parameter is derived from
the model grid. Parameters that are part of the boundary conditions of the unconfined
aquifer are listed under [Constant Head](@ref) and [Boundary conditions](@ref) of the Model
parameters section.

The wflow_sediment model consists of the vertical [Soil Erosion](@ref) concept and the input
parameters for this concept are described in the [Sediment](@ref) section of the Model
parameters. The parameters of the lateral [Sediment Flux in overland flow](@ref) concept are
described in the [Overland flow](@ref) section of the Model parameters. Parameters of this
component are not directly set by data from static input. The input parameters of the
lateral concept [River Sediment Model](@ref) are listed in [River flow](@ref) of the Model
parameters section.

The Model parameters section lists all the parameters per Model component and these Tables
can also be used to check which parameters can be part of the output, see also [Output
NetCDF section](@ref) and [Output CSV section](@ref).

For each Wflow Model a test model is available that can help to understand the data
requirements and the usage of each Model. The TOML configuration file per available model
are listed in the Table below:

|  model  | TOML configuration file |
|:--------------- | ------------------|
| [SBM + Kinematic wave](@ref) | [sbm_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_config.toml) |      
| [SBM + Groundwater flow](@ref) | [sbm\_gwf\_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_gwf_config.toml) | 
| [HBV model](@ref) | [hbv_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/hbv_config.toml) | 
| [wflow_sediment](@ref) | [sediment_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sediment_config.toml) |

The associated Model files (input static, forcing and state files) can easily be downloaded
and for this we share the following Julia code (per Model) that downloads the required files
to your current working directory:

[SBM + Kinematic wave](@ref):
```julia
# urls to TOML and NetCDF of the Moselle example model
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.4/staticmaps-moselle.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/forcing-2000.nc"
instates = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.2/instates-moselle.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "sbm_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-moselle.nc"))
download(forcing, joinpath(datadir, "forcing-moselle.nc"))
download(instates, joinpath(datadir, "instates-moselle.nc"))
download(toml_url, toml_path)
```

[SBM + Groundwater flow](@ref):
```julia
# urls to TOML and NetCDF of the Moselle example model
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.2/staticmaps-sbm-groundwater.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.1/forcing-sbm-groundwater.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "sbm_gwf_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-sbm-groundwater.nc"))
download(forcing, joinpath(datadir, "forcing-sbm-groundwater.nc"))
download(toml_url, toml_path)
```

[HBV model](@ref):
```julia
# urls to TOML and NetCDF of the Moselle example model
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.1/staticmaps-lahn.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/forcing-lahn.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "hbv_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-lahn.nc"))
download(forcing, joinpath(datadir, "forcing-lahn.nc"))
download(toml_url, toml_path)
```

[wflow_sediment](@ref):
```julia
# urls to TOML and NetCDF of the Moselle example model
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.3/staticmaps-moselle-sed.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.3/forcing-moselle-sed.nc"
instates = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/instates-moselle-sed.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "sediment_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-moselle-sed.nc"))
download(forcing, joinpath(datadir, "forcing-moselle-sed.nc"))
download(instates, joinpath(datadir, "instates-moselle-sed.nc"))
download(toml_url, toml_path)
```

For running these test model see also [Usage](@ref) and [Command Line Interface](@ref).

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