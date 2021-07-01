# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed
- Changed length units for lateral subsurface flow component from millimeter to meter. This
  means that state netCDF files from previous versions can only be reused if `ssf` is
  divided by 10^9.
- Add snow and glacier processes to wflow\_sbm figure of the documentation.

### Added
- Multi-threading of vertical SBM concept and lateral kinematic wave components (overland,
  river and subsurface flow) of wflow\_sbm model [SBM + Kinematic wave](@ref).
- Improved error message for CSV Reducer.

### Fixed
- Calculation of volumetric water content of vertical SBM (soil layers and root zone).
- Update of `satwaterdepth` in SBM (evaporation was only substracted from a local
  variable, and not from `sbm.satwaterdepth`).
- Linkedlocs index is fixes.
- Interpolation function for csv tables lakes is fixed.

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
