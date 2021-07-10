# Run from Delft-FEWS

Wflow integrates easily as part of an operational system by linking to the
[Delft-FEWS](https://oss.deltares.nl/web/delft-fews/) platform. Delft-FEWS integrates data
and models, and is for example used in many active flood forecasting systems around the
world.

This can be done without a model adapter that provides the interface between Delft-FEWS and
an external model (or module). This is possible because time information in the TOML
configuration file is optional and Delft-FEWS can import and export NetCDF files. When time
information is left out from the TOML configuration file, the `starttime`, `endtime` and
`timestepsecs` (timestep) of the run is extracted from the NetCDF forcing file by Wflow. 

To indicate that a Wflow model runs from Delft-FEWS, the following setting needs to be
specified in the main section of the TOML configuration file:

```toml
fews_run = true  # optional, default value is false
```

This ensures that Wflow offsets the time handling, to meet the expectations of Delft-FEWS.
