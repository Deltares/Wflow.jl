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
##### Dimensions #####

Name                                            Length
-------------------------------------------------------------------------
lat                                             313
time                                            366
lon                                             291

##### Variables #####

Name                        Type          Dimensions
-------------------------------------------------------------------------
lat                         FLOAT         lat
PET                         FLOAT         lon lat time
time                        INT64         time
P                           FLOAT         lon lat time
TEMP                        FLOAT         lon lat time
lon                         DOUBLE        lon

##### Attributes #####

Variable          Name              Value
-------------------------------------------------------------------------
lat               units             degrees_north
lat               long_name         latitude
lat               _FillValue        NaN
PET               _FillValue        NaN
time              units             days since 2000-01-01 00:00:00       
time              calendar          proleptic_gregorian
P                 _FillValue        NaN
TEMP              _FillValue        NaN
lon               _FillValue        NaN
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