# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed
- Model types `sbm_gwf` and `hbv` use the same approach for the calculation of the drain
  width `dw` as model type `sbm`.

### Added
- Option to use the local inertial model for river flow as part of the [SBM + Kinematic
  wave](@ref). See also [SBM + Local inertial river](@ref).
- Additional log messages and log file as output, see also [Logging](@ref).

### Fixed
- Model type `hbv`: the surface width for overland flow was not corrected with the river
  width.
- Fixed use of absolute path for `path_forcing` in TOML file, which gave an error in Wflow
  v0.5.1.
- Fixed a crash when using glaciers.
- When the surface flow width for overland flow is zero, the water level `h` of the
  kinematic wave should not be calculated, otherwise this results in `NaN` values. When the
  model is initialized from state files, `q` and `h` are set to zero for indices with a zero
  surface flow width.
- Fixed how number of iterations `its` for kinematic wave river flow are calculated during
  initialization when using a fixed sub-timestep (specified in the TOML file). For a model
  timestep smaller than the fixed sub-timestep an InexactError was thrown.
- Fixed providing a cyclic parameter when the NetCDF variable is read during model
  initialization with `ncread`, this gave an error about the size of the NetCDF `time`
  dimension.

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
  capillary rise ceases, and a coefficient `cap_n`. See also [Transpiration and soil
  evaporation](@ref). Multiplying the scaling factor with the ratio of model time step and
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
- Update of `satwaterdepth` in SBM (evaporation was only substracted from a local variable,
  and not from `sbm.satwaterdepth`).
- Fixed the index `lowerlake_ind` of `NaturalLake`. Linked lakes were not simulated
  correctly (flows between upstream and downsteam lake).
- Interpolation function `interpolate_linear(x, xp, fp)` for CSV tables lakes. When `x` was
  larger or smaller than `xp`, maximum(`xp`) or minimum(`xp`) was returned instead of
  maximum(`fp`) or minimum(`fp`).
- Fixed model timestep of reservoirs (`SimpleReservoir`) and lakes (`NaturalLake`) in
  relation to precipitation and evapotranspiration fluxes. This was set to the fixed Wflow
  `basetimestep` of 86400 s, and should be set to the actual model time step from the TOML
  configuration file.
- Add `flux` from `Drainage` (`GroundwaterFlow`) in the `sbm_gwf_model` to the overland flow
  component instead of the river component of the kinematic wave.
- Fixed option `constanthead = false` (TOML file), `constant_head` field of
  `GroundwaterFlow` was not defined in this case. Fixed this by initialising empty fields
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
- Ignore extra dimensions in input NetCDFs if they are size 1

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
  applicable for wflow\_sbm. The `f` parameter needs to be provided as part of the NetCDF
  model parameter file.
- Grid properties as cell length and elevation now stored as part of the
  `model.land.network` component and not as part of the vertical model components, as it is
  not used by these components. `altitude` (elevation) should now be provided as part of the
  `[input]` section of the TOML configuration file, and not as part of the
  `[input.vertical]` section.
- Removed parameter ``\theta_{e}`` from SBM struct (not used in update). Parameters
  ``\theta_{s}`` and ``\theta_{r}`` included separately (instead of ``\theta_{e}``) in
  `LateralSSF struct`, now directly linked to SBM parameters.
- Improve error messages (NetCDF and cyclic flow graph).

### Added
- Export of NetCDF scalar timeseries (separate NetCDF file from gridded timeseries). This
  also allows for importing these timeseries by Delft-FEWS (General Adapter).

### Fixed
- Model parameter Manning's `n` now used during the update of the `struct SurfaceFlow`,
  to change the related ``\alpha`` parameter of the kinematic wave for channel flow.
