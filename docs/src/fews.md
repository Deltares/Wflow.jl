# Run from Delft-FEWS

Wflow can be linked to the Flood forecasting system
[Delft-FEWS](https://oss.deltares.nl/web/delft-fews/), without a model adapter that provides
the interface between Delft-FEWS and an external model (or module). This is possible because
time information in the TOML configuration file is optional and Delft-FEWS can import and
export NetCDF files. When time information is left out from the TOML configuration file, the
`starttime`, `endtime` and `timestepsecs` (timestep) of the run is extracted from the NetCDF
forcing file by Wflow. 

To indicate that a Wflow model runs from Delft-FEWS, the following setting needs to be
specified in the TOML configuration file:

```toml
fews_run = true                         # optional, default value is false
```