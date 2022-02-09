# Step 1: Understanding the requirements

In order to run wflow, several files are required. These consist of a settings file and input data. The input data is typically separated into static maps and forcing data, both are supplied in a NetCDF file. A brief overview of the different files:

 - The `settings.toml` file contains information of the simulation period, links to the input files (and their names in the NetCDF files), and links the correct names of the variables in the NetCDF files to the variables and parameters of wflow.
 - The `staticmaps.nc` file contains spatial information on the elevation, locations of the gauges, land-use, drainage direction, etc. This can also contain maps with parameter values.
 - The `forcing.nc` file contains the precipitation, temperature and potential evaporation time series (as a 3D array). 

There are several model configurations supported by wflow. These model configurations required slightly different input requirements, yet the general structure is similar for each model. In the following pages, some examples will be given on how to prepare a basic `wflow\_sbm` model. Example data for other model configurations is provided in the section with [sample data](@ref sample_data). The following model configurations are supported in wflow:

 - wflow\_sbm + kinematic wave
 - wflow\_sbm + groundwater flow
 - wflow\_hbv + kinematic wave
 - wflow\_sediment as post processing of wflow_sbm or wflow_hbv output

