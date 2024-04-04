# Step 3: Preparing the input data

As mentioned before, the input data can be classified into two types:

 - Meteorological forcing: maps with timeseries for each model pixel, with values for
   precipitation, temperature, and potential evaporation. This data should be provided as a
   three-dimensional dataset, with the `x`, `y` and `time` dimensions.
 - Static maps.


## Meteorological data

Meteorological data is provided as a single netCDF file, with several variables containing
the forcing data for precipitation, temperature and potential evaporation. The code snippet
below shows the contents of the example file (downloaded [here](@ref wflow_sbm_data)), and
displaying the content with `NCDatasets` in Julia. As can be seen, each forcing variable
(`precip`, `pet` and `temp`) consists of a three-dimensional dataset (`x`, `y`, and `time`),
and each timestep consists of a two-dimensional map with values at each gridcell. Only
values within the basin are required.

```
Group: /

Dimensions
   time = 366
   y = 313
   x = 291

Variables
  time   (366)
    Datatype:    Int64
    Dimensions:  time
    Attributes:
     units                = days since 2000-01-02 00:00:00
     calendar             = proleptic_gregorian

  y   (313)
    Datatype:    Float64
    Dimensions:  y
    Attributes:
     _FillValue           = NaN

  x   (291)
    Datatype:    Float64
    Dimensions:  x
    Attributes:
     _FillValue           = NaN

  spatial_ref
    Attributes:
     crs_wkt              = GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
     x_dim                = x
     y_dim                = y
     dim0                 = time

  precip   (291 × 313 × 366)
    Datatype:    Float32
    Dimensions:  x × y × time
    Attributes:
     _FillValue           = NaN
     unit                 = mm
     precip_fn            = era5
     coordinates          = idx_out spatial_ref mask

  idx_out   (291 × 313)
    Datatype:    Int32
    Dimensions:  x × y

  mask   (291 × 313)
    Datatype:    UInt8
    Dimensions:  x × y

  pet   (291 × 313 × 366)
    Datatype:    Float32
    Dimensions:  x × y × time
    Attributes:
     _FillValue           = NaN
     unit                 = mm
     pet_fn               = era5
     pet_method           = debruin
     coordinates          = idx_out spatial_ref mask

  temp   (291 × 313 × 366)
    Datatype:    Float32
    Dimensions:  x × y × time
    Attributes:
     _FillValue           = NaN
     unit                 = degree C.
     temp_fn              = era5
     temp_correction      = True
     coordinates          = idx_out spatial_ref mask

Global attributes
  unit                 = mm
  precip_fn            = era5
```

!!! note
    Wflow expects right labeling of the forcing time interval, e.g. daily precipitation
    at 01-02-2000 00:00:00 is the accumulated total precipitation between 01-01-2000
    00:00:00 and 01-02-2000 00:00:00.


## Static data


### List of essential static data

The list below contains a brief overview of several essential static maps required to run
wflow. These NC variables names refer to the example data of the wflow\_sbm + kinematic wave
model (see [here](@ref wflow_sbm_data)). Example data for the other model configurations can
be found [here](@ref sample_data).

Description | NC variable name | unit
--- | --- | ---
Flow direction (1-9) | `wflow_ldd` | -
Map indicating the river cells (0-1) | `wflow_river` | -
The length of the river | `wflow_riverlength` | m
The width of the river | `wflow_riverwidth` | m
Mask of the basin | `wflow_subcatch` | -
Land slope | `Slope` | m m$^{-1}$
River slope | `RiverSlope` | m m$^{-1}$

As mentioned before, the model parameters can also be defined as spatial maps. They can be
included in the same netCDF file, as long as their variable names are correctly mapped in
the TOML settings file. See the section on [example models](@ref sample_data) on how to
use this functionality.
