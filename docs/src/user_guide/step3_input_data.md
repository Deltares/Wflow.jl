# Step 3: Preparing the input data

As mentioned before, the input data can be classified into two types:

 - Meteorological forcing: maps with timeseries for each model pixel, with values for
   precipitation, temperature, and potential evaporation. This data should be provided as a
   three-dimensional dataset, with the `x`, `y` and `time` dimensions. 
 - Static maps.


## Meterological data

Meteorological data is provided as a single NetCDF file, with several variables containing
the forcing data for precipitation, temperature and potential evaporation. The code snippet
below shows the contents of the example file (downloaded [here](@ref wflow_sbm_data)), and
opened with `xarray` in Python. As can be seen, each variable consists of a
three-dimensional dataset (`lon`, `lat`, and `time`), and each timestep consists of a
two-dimensional map with values at each gridcell. Only values within the basin are required.

```
<xarray.Dataset>
Dimensions:  (lon: 291, lat: 313, time: 366)
Coordinates:
  * lon      (lon) float64 5.429 5.438 5.446 5.454 ... 7.821 7.829 7.837 7.846
  * lat      (lat) float32 50.42 50.41 50.4 50.4 ... 47.85 47.84 47.83 47.82
  * time     (time) datetime64[ns] 2000-01-01 2000-01-02 ... 2000-12-31
Data variables:
    P        (time, lat, lon) float32 ...
    TEMP     (time, lat, lon) float32 ...
    PET      (time, lat, lon) float32 ...
```

## Static data


### List of essential static data

The list below contains a brief overview of several essential static maps required to run
wflow. These NC variables names refer to the example data of the wflow\_sbm model (see
[here](@ref wflow_sbm_data)).

Description | NC variable name | unit
--- | --- | ---
Altitude | `wflow_dem` | m
Flow direction (1-9) | `wflow_ldd` | -
Location of cells that are pits | `wflow_pits` | -
Map indicating the river cells (0-1) | `wflow_river` | -
The length of the river | `wflow_riverlength` | m
The width of the river | `wflow_riverwidth` | m
Mask of the basin | `wflow_subcatch` | -
Land slope | `Slope` | m m$^-1$
River slope | `RiverSlope` | m m$^-1$

As mentioned before, the model parameters can also be defined as spatial maps. They can be
included in the same NetCDF file, as long as their variable names are correctly mapped in
the TOML settings file.