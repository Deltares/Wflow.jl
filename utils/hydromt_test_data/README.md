# Rebuilding the Wflow test data with HydroMT

This folder holds a pipeline that rebuilds the test models of `Wflow/test`, i.e. the TOML
files and the files in `Wflow/test/data/input`, with
[HydroMT-Wflow](https://deltares.github.io/hydromt_wflow/latest/). It documents how the
test data is constructed and makes it possible to regenerate it, for example after
changes in the input data or in HydroMT-Wflow. The existing test data
(`pixi run download-test-data`) is not changed by the pipeline.

The generated data is not identical to the existing test data, which was built with older
versions of HydroMT-Wflow, data sources and Wflow. Replacing the test data therefore
requires updating the reference values of the tests.

## Requirements

- Access to the Deltares data catalog (`deltares_data`, P: drive). Other catalogs can be
  used with `--data-libs`, as long as they provide the data sources of the workflows.
- Wind speed and net radiation (land surface temperature forcing) are derived from hourly
  ERA5 in `era5.zarr` on the P drive (`data_catalog.yml`, chunked for regional time
  series). Net radiation follows FAO-56: net shortwave radiation with an albedo of 0.23
  and estimated net longwave radiation, because ERA5 net thermal radiation is not
  available on the P drive (see `era5_lst_forcing` in `steps.py`). ERA5 net solar
  radiation can be used with `ssr_source: era5_hourly`, but those files are unchunked
  global arrays, which makes reading a regional time series slow. To use ERA5 net
  radiation (ssr + str), set
  `source: earthdatahub` ([DestinE Earth Data Hub](https://earthdatahub.destine.eu),
  requires a free API key in `~/.netrc`) or `source: arco` (Google ARCO-ERA5, no key,
  but chunked per hour for the whole globe, ~70 GB for a year, so limit the period with
  `start` and `end`) for `pipeline.era5_lst_forcing` in `workflows/moselle_sbm.yml`.
- Internet access for the HydroMT-Wflow water demand example data (Piave region).
- The `hydromt` pixi environment (Python, HydroMT and HydroMT-Wflow) and the default
  pixi environment (Julia) for the Wflow runs.

## Usage

From the repository root:

```sh
pixi run -e hydromt build-test-data                 # all stages
pixi run -e hydromt build-test-data build export    # selected stages
pixi run -e hydromt build-test-data --only piave_sbm piave_spinup
```

The output is written to `utils/hydromt_test_data/output` (or `--output DIR`):

- `models/`: the HydroMT models with HydroMT-Wflow names.
- `staging/`: the TOML files of the HydroMT models after renaming (see below).
- `test/`: the generated test files, with the same layout as `Wflow/test`.
- `runs/`: the TOML files and output of the Wflow runs.

## Stages

1. `build`: HydroMT-Wflow builds the models in `workflows/*.yml` (`moselle_sbm`,
   `moselle_sediment` as an update of `moselle_sbm`, `piave_sbm` and `rhine_sbm`). Steps
   prefixed with `pipeline.` are defined in `steps.py` and add data that HydroMT-Wflow
   does not provide, e.g. constant parameter maps, land surface temperature parameters and
   forcing, the maps of the `sbm_gwf` model type and the test-specific maps used in
   `Wflow/test/*.jl`.
2. `export`: renames the HydroMT-Wflow variables to the netCDF variable names of the tests
   (`naming/*.yml`) and writes the static maps and forcing (`export` section of the
   workflows). It also writes the data that is not built with HydroMT:
   - the synthetic `sbm_gwf` test case (`synthetic.py`, `synthetic_sbm_gwf_forcing.csv`),
   - `forcing-calendar-noleap.nc`, a subset of the Moselle forcing with a noleap calendar,
   - the reservoir storage and rating curve tables (`reservoir_curves/`).
3. `configs`: writes the test TOML files from `configs/*.yml` and the Wflow run TOML files
   from `runs/*.yml`. The `inputs` of a specification are taken from the exported HydroMT
   model config (`model`), all other TOML content is given in `config` (dotted keys).
   Run specifications can `extend` a test specification and `remove` keys.
4. `run`: runs Wflow to produce the initial states (spin-up runs with a cold start), the
   forcing of the sediment model (gridded output of the Moselle sbm model) and the
   overland flow kinematic wave test data of the Rhine (`staticmaps-rhine.nc`).
5. `compare`: reports the differences between the generated files and `Wflow/test`.

## Overview of the test data

| File | Source |
| --- | --- |
| `staticmaps-moselle.nc`, `forcing-moselle.nc` | `workflows/moselle_sbm.yml` |
| `staticmaps-moselle-sed.nc` | `workflows/moselle_sediment.yml` |
| `forcing-moselle-sed.nc` | `runs/moselle_sediment_forcing.yml` |
| `staticmaps-piave.nc`, `forcing-piave.nc` | `workflows/piave_sbm.yml` |
| `staticmaps-rhine.nc` | `workflows/rhine_sbm.yml`, `runs/rhine_kinematic_wave.yml` |
| `instates-*.nc` | `runs/*_spinup.yml` |
| `forcing-calendar-noleap.nc` | subset of `forcing-moselle.nc` |
| `staticmaps-sbm-groundwater.nc`, `forcing-sbm-groundwater-part*.nc` | `synthetic.py` |
| `reservoir_*.csv` | `reservoir_curves/` |

## Known differences with the existing test data

- The time stamp of `instates-moselle.nc` is the end of a one year spin-up (2000-12-31)
  instead of 2000-01-02.
- `h_av_land` and `h_av_river` in the states are copies of `h_land` and `h_river`.
- `CanopyGapFraction` (`sediment_config.toml`) is not in the static maps, as in the
  existing test data.
- The land surface temperature parameters are only defined for the GlobCover land use
  classes of the Moselle model (`lst_mapping.csv`).
- Net radiation in `forcing-moselle.nc` is the FAO-56 net radiation instead of ERA5 net
  solar plus net thermal radiation (`ssr + str`), see above.
- Constant parameters are written as maps (`pipeline.constant_maps`), because HydroMT-Wflow
  v1 writes them as values in the TOML file, while the tests read them from the static maps.
- The reservoir areas in the current Deltares reservoir database are smaller than in the
  existing test data (Piave: 0.59 and 0.62 km2 instead of 1.79 and 1.31 km2), which
  affects the derived reservoir parameters.
- The Moselle sbm and Piave models were built with an older HydroMT version, whose model
  grid is not aligned to multiples of the model resolution. `pipeline.basemaps` reproduces
  these grids (`grid_bounds`): the hydrography is clipped to these bounds extended by whole
  model cells, and border rows or columns without active cells are padded.
  The sediment model is an update of the Moselle sbm model
  and shares its grid, while the existing sediment test data is on the aligned grid
  (shifted by 3 and 4 hydrography cells in x and y).
- The Rhine model (`staticmaps-rhine.nc`) covers a smaller domain (774 x 678 instead of
  927 x 700 cells); the existing file has no coordinates to reproduce its extent.
- The Piave land use (GLCNMO from the HydroMT-Wflow water demand example data) differs
  from the existing test data on about 2000 cells: the existing data has classes (14, 16)
  that are (almost) absent in the example data. This changes the paddy and non-paddy
  irrigation areas and the land use based subcatchments of the Piave masking test.
