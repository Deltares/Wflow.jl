# Step 1: Understanding the requirements

To run wflow, several files are required. These include a settings file and input data. The
input data is typically separated into static maps and forcing data, and both are provided
in netCDF files, except for lake storage and rating curves that are supplied via CSV files.
Below is a brief overview of the different files:

 - The `settings.toml` file contains information on the simulation period, links to the
   input files (and their names in the netCDF files), and connect the correct variable names
   in the netCDF files to the variables and parameters of wflow.
 - The `staticmaps.nc` file contains spatial information such as elevation, guage locations,
   land use, and drainage direction, etc. This file can also contain maps with parameter
   values.
 - The `forcing.nc` file contains time series data for  precipitation, temperature and
   potential evaporation (as a 3D array).

Wflow supports several model configurations, each requiring slightly different input, but
with a similar general structure. A wflow model configuration consists of a `vertical`
concept like [SBM](@ref vert_sbm), [HBV](@ref vert_hbv) or [FLEXTOPO](@ref vert_flextopo) in
combination with `lateral` concepts that control how water is routed for example over the
land or river domain. For the wflow\_sbm model, different model configurations are possible.
The following configurations are supported in wflow:

 - wflow\_sbm:
    - SBM + kinematic wave for subsurface and surface flow
    - SBM + kinematic wave for subsurface and overland flow + local inertial river (+
      optional floodplain)
    - SBM + kinematic wave for subsurface flow + local inertial river (1D) and land (2D)
    - SBM + groundwater flow + kinematic wave for surface flow
 - wflow\_hbv: HBV + kinematic wave for surface routing
 - wflow\_flextopo: FLEXTOPO + kinematic wave for surface routing
 - wflow\_sediment as post processing of wflow\_sbm or wflow\_hbv output

In the following pages, some examples will be given on how to prepare a basic wflow\_sbm
model. Sample data for other model configurations is provided in the [sample data](@ref
sample_data) section.