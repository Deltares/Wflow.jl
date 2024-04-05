# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Fixed
- Added missing BMI function `get_grid_size`, it is used for unstructured grids, for example
  to get the length of arrays returned by BMI functions `get_grid_x` and `get_grid_y`.
- Added a check for the solution of the quadratic equation as part of the Modified Puls
  approach for lake outflow. Lower limit should be zero (very small negative values can
  occur).
- Limit lake evaporation (added variable `actevap`) and lake outflow to prevent negative
  lake storage. The variable `actevap` has also been added to the reservoir module.
- The `set_states` function for model type `sbm` with local inertial routing for river and
  land component.
- Inflow to reservoir and lake locations for local inertial routing: 1) with floodplain
  routing, the floodplain discharge was not added to the inflow of these locations, 2) the
  `to_river` variable from overland flow and lateral subsurface flow was not added to the
  inflow of these locations.
- Close netCDF `NCDataset` with state variables in extended BMI function `save_state`.
- `BMI.update_until` could throw an `InexactError: Int64` caused by a not whole number. This
  is fixed by applying `round()`.

### Changed
- Stop exposing scalar variables through BMI. The `BMI.get_value_ptr` function was
  not working correctly for scalar model variables (a `view` was applied). Only a few scalar
  model parameters are defined, and it is not expected that exposing these variables is
  required (e.g. for model coupling) while code changes for these variables (including
  struct fields) are required.
- The local inertial routing (constant) river boundary condition `depth` at a ghost node
  (downstream river outlet) was set at 0.0 in the code and can now be provided through the
  model parameter netCDF file (`riverdepth_bc`). In addition, the boundary condition river
  length `dl` at a ghost node should be set through the model parameter netCDF file
  (`riverlength_bc`), providing this boundary condition through the `[model]` settings in
  the TOML file (`model.riverlength_bc`) has been deprecated.

### Added
- Total water storage as an export variable for `SBM` concept. This is the total water stored
  per grid cell in millimeters. Excluded from this variable are the floodplain, lakes and
  reservoirs.
- Checks to see if all states are covered in the .toml file. If not all states are covered,
  an error is thrown. If there are more states specified than required, these states are
  ignored (with a warning in the logging), and the simulation will continue.
- Support different vertical hydraulic conductivity profiles for the `SBM` concept. See also
  the following sections: [The SBM soil water accounting scheme](@ref) and [Subsurface flow
  routing](@ref) for a short description.
- Local inertial routing to `sbm_gwf` model type.

## v0.7.3 - 2024-01-12

### Fixed
- Documentation: add leakage term to the wflow\_sbm figure, document external input
  parameter `ksathorfrac` and fix description of adding external `inflow` to the kinematic
  wave.
- Fixed BMI functions (e.g. `BMI.get_value`) that deviated from BMI specifications
  (BasicModelInterface.jl), including function arguments, return types and the BMI
  specification that arrays are always flattened (this was not the case for variables stored
  as 2-dimensional arrays or as vector of SVectors).
- Bump compat for NCDatasets to 0.13, 0.14.
- The solution for lake outflow as part of the Modified Puls Approach. The inflow and
  outflow variables are defined for period `Δt`, and not at `t1` and `t2` (instantaneous) as
  in the original mass balance equation of the Modified Puls Approach. Because of this, the
  terms of the quadratic equation (and solution) were fixed.
- Use `kvfrac` for the computation of vertical saturated hydraulic conductivity at the
  bottom of the soil layer, since `kvfrac` is also used for the computation of vertical
  unsaturated flow.

### Changed
- For cyclic parameters different cyclic time inputs are supported (only one common cyclic
  time (for example daily or monthly) was allowed).
- BMI: 1) added grid information (type and location) and whether a variable can be exchanged
  to metadata Structs, 2) extend model grid functions for wflow components that store
  variables on `edges` (local inertial model) with `get_grid_edge_count` and
  `get_grid_edge_nodes`.

### Added
- Functions for loading and saving states and getting the `starttime` in Unix time format.
  This is convenient for coupling with OpenDA (and other external) software. The set states
  functionality from the initialization function has been moved to a separate `set_states`
  function for each `Model` type, to support loading states independent of initialization.

## v0.7.2 - 2023-09-27

### Fixed
- Water balance of modified Rutter interception model. The sum of the stemflow partitioning
  coefficient `pt` and free throughfall coefficient `p` could get larger than 1, resulting
  in an overestimation of stemflow and throughfall amounts and negative net interception
  amounts. And the first drainage amount `dd`, controlled by a change over time in canopy
  storage capacity `cmax`, should not be subtracted from precipitation to compute net
  interception.
- The `netinterception` (net interception) computed by the modified Rutter interception
  model was stored as `interception` in `SBM`, while this should be the `interception`
  (interception loss by evaporation) output of the modified Rutter interception model. The
  `interception` of `SBM` is used to compute the total actual evapotranspiration `actevap`.

## v0.7.1 - 2023-06-30

### Fixed
- State initialization of 1D floodplain `volume`. In the initialization function the wrong
  field name of type `FloodPlainProfile` was used (`area` instead of `a`).

## v0.7.0 - 2023-06-12

### Fixed
- `BMI.get_time_units` now gets called on the model rather than the type, like all other BMI
  functions, except `BMI.initialize`. Also it returns "s" instead of "seconds since
  1970-01-01T00:00:00", in line with the BMI specification.
- Added the `interception` component to total actual evapotranspiration `actevap` of `SBM`
  (was defined as the sum of soil evaporation, transpiration and open water evaporation).

### Changed
- The time values returned in the BMI interface are no longer in seconds since 1970, but in
  seconds since the model start time. This is more in line with standard BMI practices.
- The `starttime` was defined one model timestep `Δt` ahead of the actual model time (the
  initial conditions timestamp (state time)). As a consequence this was also the case for
  the current model time. To allow for an easier interpretation of wflow time handling,
  either through BMI or directly, the `starttime` is now equal to the state time, resulting
  in current model times without an offset.
- Using more than 8 threads can result in too much overhead with `Threads.@threads`. After
  performance testing, this has been changed for kinematic wave routing and the vertical
  `SBM` concept to spawning tasks with `Threads@spawn` for number of threads <= 8, where
  each task iterates over a chunk of size `basesize`. For more than 8 threads the low
  overhead threading `Polyester.@batch` (including the `minbatch` argument) is used. For
  local inertial routing the use of `Threads.@threads` has been changed to threaded loop
  vectorization (river and 1D floodplain local inertial momentum equation) and
  `Polyester.@batch`.

### Added
- For (regulated) lakes with rating curve of type 1 (H-Q table), lake `storage` above the
  `maximumstorage` (based on maximum water level from the H-Q table) is spilled
  instantaneously (overflow) from the lake.
- Added support to use `sum` as a reducer function for csv and scalar output options.

## v0.6.3 - 2023-03-01

### Fixed
- Removed error when `_FillValue` is present in the time dimension of the forcing netCDF
  file. The simulation is allowed to continue with the attribute present, given that there
  are no missing values in the time dimension. This is checked by the code, and an error is
  thrown if this is the case.
- Column index of daily lake rating curves. This was incorrectly based on `dayofyear` with a
  maximum of 365. The column index should be based on julian day (leap days are not
  counted).

### Changed
- `NCDatasets` version. Reading the `time` dimension of multifile netCDF file became very
  slow since `NCDatasets` v0.12.4, this issue has been solved in v0.12.11.
- Store the `time` dimension of the forcing netCDF file as part of the struct `NCreader`
  instead of calling `dataset["time"][:]` each time step when loading forcing data.

### Added
- Show total duration of simulation in the log file (info), and show the current time at
  execution of each timestep (debug).
- Support for exponential decline in horizontal conductivity in the sbm\_gwf concept. This can
  be enabled using the `conductivity_profile` setting (either "uniform" or "exponential"). If
  set to "exponential", it exponentially reduces the `kh0` (or `conductivity`) based on the
  value of `gwf_f` to the actual horizontal conductivity (`k`).
- An optional 1D floodplain schematization for the river flow inertial model, based on
  provided flood volumes as a function of flood depth per river cell. See also the following
  sections: [SBM + Local inertial river and floodplain](@ref) and [River and floodplain
  routing](@ref) for a short description, and the following section for associated [model
  parameters](@ref local-inertial_floodplain_params).

## v0.6.2 - 2022-09-01

### Fixed
- Two issues related to reservoir and lake locations as part of local inertial model: 1)
  added as boundary points to the update of overland flow, 2) fixed check reservoir and lake
  location in update river flow.
- Limit flow (`qx`, `qy` and `q`) in local inertial model when water is not available (set
  to zero).
- In the grid output netCDFs, don't set the `_FillValue` attribute, since the CF conventions
  don't allow it.
- Glacier parameter `g_sifrac` needs to be converted during initialization (time dependent).
- The check that the sum of adaptive timesteps (`Δt`) of the local inertial model (1D and
  2D) does not exceed the model timestep.

### Changed
- Changed depth `h` for reservoir and lake locations as part of the river local inertial
  model from `bankfull_depth` to zero.

### Added
- External inflow to the [SBM + Local inertial river (1D) and land (2D)](@ref) model
  configuration.

## v0.6.1 - 2022-04-26

### Fixed
- Fixed an error with the log file, when writing to a folder that does not (yet) exists.
  Now, the folder is created prior to writing the log file.
- Fixed a MethodError for `read_dims`, thrown when reading netCDF data with NCDatasets.jl
  0.12.3 or higher.

## v0.6.0 - 2022-04-14

### Added
- The [FLEXTopo](@ref config_flextopo) model.
- Set a (different) uniform value for each index of input parameters with an extra
  dimension. A list of values (instead of one uniform value) that should be equal to the
  length of the extra dimension can be provided in the TOML file. See also [Modify
  parameters](@ref).
- Modify input parameters with an extra dimension at multiple indices with different `scale`
  and `offset` values, through the TOML file. See also [Modify parameters](@ref).
- Added support for netCDF compression for gridded model output, through the option
  `compressionlevel` in the `[output]` section in the TOML file. The setting defaults to 0
  (no compression), but all levels between 0 and 9 can be used.

### Changed
- Re-organized the documentation: moved explanation of different model concepts to a model
  documentation section, added a user guide to explain setting up the model, added new
  figures to the description of wflow\_sbm.
- The unit of lateral subsurface flow `ssf` of `LateralSSF` is now ``m^3 d^{-1}``. The unit
  was ``m^3 t^{-1}``, where ``t`` is the model timestep. Other flow variables are already
  stored independently from ``t``, this allows for easier interpretation and to use states
  independently of ``t``.
- Changed the reference level of water depth `h` and `h_av` of 2D overland flow
  (`ShallowWaterLand`) for cells containing a river from river bed elevation `zb` to cell
  elevation `z`.

### Fixed
- Fixed calculation of average water depth `h_av` of 2D overland flow (`ShallowWaterLand`)
  with the local inertial approach. The summation of `h` was not correct, resulting in too
  low values for `h_av`. For river cells of 2D overland flow `h_av` was only updated as part
  of the river domain (`ShallowWaterRiver`), this value is now also updated as part of the
  land domain (`ShallowWaterLand`).
- Fixed the following two flow width issues for 2D overland flow of the local inertial
  model: 1) The flow widths `xwidth` and `ywidth` were mapped incorrectly (reversed) to the
  flow calculation in `x` and `y` direction, respectively. 2) For river cells the effective
  flow width for overland flow was not determined correctly: the `riverwidth` vector
  supplied to the function `set_effective_flowwidth!` covered the complete model domain and
  should have covered the river domain, resulting in incorrect river widths and thus
  incorrect effective flow widths for river cells.

## v0.5.2 - 2022-02-03

### Changed
- Model types `sbm_gwf` and `hbv` use the same approach for the calculation of the drain
  width `dw` as model type `sbm`.
- Renamed `h_bankfull` parameter to `bankfull_depth` for consistency between kinematic-wave
  and local-inertial river routing. When using the old name under the
  `[input.lateral.river]` section of the TOML file, it will work but it is suggested to
  update the name.

### Added
- Additional log messages and log file as output, see also [Logging](@ref logging_toml).
- Option to use the local inertial model for river flow as part of the `sbm` model type. See
  also [SBM + Local inertial river](@ref).
- Option to use the local inertial model for 1D river flow combined with 2D overland flow as
  part of the `sbm` model type. See also [SBM + Local inertial river (1D) and land
  (2D)](@ref).

### Fixed
- Model type `hbv`: the surface width for overland flow was not corrected with the river
  width.
- Fixed use of absolute path for `path_forcing` in TOML file, which gave an error in wflow
  v0.5.1.
- Fixed a crash when glacier processes are simulated as part of the `hbv` concept (Δt was
  not defined).
- When the surface flow width for overland flow is zero, the water level `h` of the
  kinematic wave should not be calculated, otherwise this results in `NaN` values. When the
  model is initialized from state files, `q` and `h` are set to zero for indices with a zero
  surface flow width.
- Fixed how number of iterations `its` for kinematic wave river flow are calculated during
  initialization when using a fixed sub-timestep (specified in the TOML file). For a model
  timestep smaller than the fixed sub-timestep an InexactError was thrown.
- Fixed providing a cyclic parameter when the netCDF variable is read during model
  initialization with `ncread`, this gave an error about the size of the netCDF `time`
  dimension.
- Fixed CSV and netCDF scalar output of variables with dimension `layer` (`SVector`).

## v0.5.1 - 2021-11-24

### Fixed
- Fixed calculation of `exfiltwater` as part of the `sbm_gwf` model type. This was based
  directly on groundwater head above the surface level, without multiplying by the
  `specific_yield`, resulting in an overestimation of `exfiltwater`. This is required since
  the groundwater model estimates the head by dividing the volume by the specific yield or
  storativity of the aquifer. So, should the groundwater table rise above surface level, the
  head above surface level does not represent a water column one to one. (This also means
  the groundwater model (slightly) overestimates heads when the head rises above the surface
  level. However, this water is immediately removed, and the head will be set to surface
  level.)

### Added
- Optional `dir_input` and `dir_output` keys in the TOML, which can be used to quickly
  change the path for all input or output files that are given as a relative path.

## v0.5.0 - 2021-11-12

### Changed
- Scaling of potential capillary rise is replaced by a common approach found in literature,
  based on the water table depth `zi`, a maximum water depth `cap_hmax` beyond which
  capillary rise ceases, and a coefficient `cap_n`. See also [Capillary rise](@ref).
  Multiplying the scaling factor with the ratio of model time step and
  `basetimestep` in the original approach resulted in (much) smaller capillary fluxes at
  sub-daily model time steps compared to daily model time steps, and is not used in the new
  approach. Parameters `cap_hmax` and `cap_n` can be set through the TOML file, parameter
  `capscale` of the previous approach is not used anymore.

### Fixed
- Conversion of `GroundwaterFlow` boundaries [``m^3 d^{-1}``] as part of model concept
  `sbm_gwf` to ``m^3 s^{-1}`` for sub-daily model time steps. For the conversion the
  `basetimestep` (86400 s) should be used (and not the model time step).

## v0.4.1 - 2021-11-04

### Changed
- The ``\alpha`` parameter of the kinematic wave has a fixed value now and is not updated
  because of changes in water height (this could result in large water balance errors). See
  also [Surface routing](@ref).
- Cyclic input for other structs than vertical are also now supported (for example cyclic
  inflow to the river).
- Moved `update_forcing!` and `update_cyclic!` functions to the `run` function. It is now
  easier to implement a custom `run` function with custom loading of input data (forcing and
  cyclic parameters).

### Added
- Check if reservoirs and lakes have downstream nodes. Without downstream nodes is not
  supported and in that case an error message is thrown that is easier to understand than
  the previous one: "ArgumentError: Collection is empty, must contain exactly 1 element."
- For the `input.path_forcing` TOML setting we use
  [Glob.jl](https://github.com/vtjnash/Glob.jl/) to expand strings like
  `"data/forcing-year-*.nc"` to a set of netCDF files that are split in time.
- External water inflow (supply/abstractions) added to the kinematic wave `inflow` in m3/s.
  It is added/removed to `sf.qlat[v]` before computing the new `q[v]` with the kinematic
  wave equation.
- Fixed values for forcing parameters are supported, see also [Fixed forcing values](@ref).

### Added
- Option to use the local inertial model for river flow as part of the [SBM + Kinematic
  wave](@ref). See also [SBM + Local inertial river](@ref).

### Fixed
- River inflow for reservoirs and lakes in the kinematic wave. This inflow was based on
  `sf.q[v]` at the previous time step, and this has been fixed to the current time step.

## v0.4.0 - 2021-09-02

### Changed
- Changed length units for lateral subsurface flow component from millimeter to meter. This
  means that state netCDF files from previous versions can only be reused if `ssf` is
  divided by 10^9.
- Add snow and glacier processes to wflow\_sbm figure of the documentation.

### Added
- Multi-threading of vertical SBM concept and lateral kinematic wave components (overland,
  river and subsurface flow) of wflow\_sbm model [SBM + Kinematic wave](@ref).
- Improved error message for CSV Reducer.
- The TOML keys `kv₀`, `θᵣ` and `θₛ` have been replaced with the ASCII versions `kv_0`,
  `theta_r` and `theta_s`, to avoid encoding issues with certain text editors. The old keys still
  work as well.

### Fixed
- Calculation of volumetric water content of vertical SBM (soil layers and root zone).
- Update of `satwaterdepth` in SBM (evaporation was only subtracted from a local variable,
  and not from `sbm.satwaterdepth`).
- Fixed the index `lowerlake_ind` of `NaturalLake`. Linked lakes were not simulated
  correctly (flows between upstream and downstream lake).
- Interpolation function `interpolate_linear(x, xp, fp)` for CSV tables lakes. When `x` was
  larger or smaller than `xp`, maximum(`xp`) or minimum(`xp`) was returned instead of
  maximum(`fp`) or minimum(`fp`).
- Fixed model timestep of reservoirs (`SimpleReservoir`) and lakes (`NaturalLake`) in
  relation to precipitation and evapotranspiration fluxes. This was set to the fixed wflow
  `basetimestep` of 86400 s, and should be set to the actual model time step from the TOML
  configuration file.
- Add `flux` from `Drainage` (`GroundwaterFlow`) in the `sbm_gwf_model` to the overland flow
  component instead of the river component of the kinematic wave.
- Fixed option `constanthead = false` (TOML file), `constant_head` field of
  `GroundwaterFlow` was not defined in this case. Fixed this by initializing empty fields
  (Vector) for struct `ConstantHead`.
- Fixed return max(0, boundary.flux[index]) to return max(0, flux) the flux should be returned
  when cell is dry, no negative value allowed.
- Fixed path to initialize lake to: dirname(static_path)
- Fixed outflow = 0, when lake level is below lake threshold. Before a negative value
  could enter the log function and model would fail.
- Fixed the lake storage initialization. For continuation runs (`reinit = false`), this
  caused the lake to be reset to the initial conditions instead of the saved state.

## v0.3.1 - 2021-05-19

### Fixed
- Ignore extra dimensions in input netCDFs if they are size 1

## v0.3.0 - 2021-05-10

### Changed
- New deposition process for coarse sediment in the reservoirs with a new parameter
  `restrapefficiency` in the sediment model.
- New variables added to the `LandSediment` and `RiverSediment` structs in order to save
  more output from the sediment model.
- Added variables `volume` and `inwater` to `SurfaceFlow` struct, this is convenient for the
  coupling with the water quality model Delwaq.
- River water level (`h`) and discharge (`q`) forced directly into the `RiverSediment`
  struct (instead of using the `OverlandFlowSediment` struct first).
- Require Julia 1.6 or later.

### Added
- Modify model parameters and forcing through the TOML file (see [Modify parameters](@ref)).
- Run wflow\_sbm (SBM + kinematic wave) in two parts (until recharge and after subsurface
  flow) from BMI, including the option to switch off the lateral subsurface component of
  wflow\_sbm.
- Support more netCDF dimension and axis order variants.

### Fixed
- Corrected a bug in sediment deposition in the river (case when incoming sediment load is
  more than the river transport capacity).
- Fixed update of `snow` and `glacierstore` fields of HBV and SBM concepts by the
  `glacier_hbv` module.

## v0.2.0 - 2021-03-26

### Changed
- Removed dependency of the `f` model parameter of wflow\_sbm on the parameters
  ``\theta_{s}``, ``\theta_{r}`` and ``M``. This approach is used in Topog\_SBM, but not
  applicable for wflow\_sbm. The `f` parameter needs to be provided as part of the netCDF
  model parameter file.
- Grid properties as cell length and elevation now stored as part of the
  `model.land.network` component and not as part of the vertical model components, as it is
  not used by these components. `altitude` (elevation) should now be provided as part of the
  `[input]` section of the TOML configuration file, and not as part of the
  `[input.vertical]` section.
- Removed parameter ``\theta_{e}`` from SBM struct (not used in update). Parameters
  ``\theta_{s}`` and ``\theta_{r}`` included separately (instead of ``\theta_{e}``) in
  `LateralSSF struct`, now directly linked to SBM parameters.
- Improve error messages (netCDF and cyclic flow graph).

### Added
- Export of netCDF scalar timeseries (separate netCDF file from gridded timeseries). This
  also allows for importing these timeseries by Delft-FEWS (General Adapter).

### Fixed
- Model parameter Manning's `n` now used during the update of the `struct SurfaceFlow`,
  to change the related ``\alpha`` parameter of the kinematic wave for channel flow.
