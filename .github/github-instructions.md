# Copilot instructions

This document provides essential information for AI coding assistants working on the Wflow project.

## Project overview

Wflow is a Julia package that provides a hydrological modeling framework, as well as several different vertical and lateral concepts that can be used to run hydrological simulations. 

- **Primary Language**: Julia (core simulation engine)
- **Domain**: Water resources modeling, hydrology, scientific computing

## Repository Structure

```
├── Wflow/                  # Core Julia package
│   ├── src/                # Source code (models, routing, config, BMI, I/O)
│   ├── test/               # Julia unit tests
│   └── Project.toml        # Julia package configuration
├── build/                  # Build tooling
│   ├── create_binaries/    # Scripts to compile standalone binaries
│   └── wflow_cli/          # CLI application package
├── server/                 # ZMQ-based BMI server
│   ├── src/                # Server source code
│   └── test/               # Server tests
├── docs/                   # Documentation (Quarto-based)
│   ├── model_docs/         # Model concept documentation
│   ├── user_guide/         # User guide pages
│   ├── getting_started/    # Installation and quickstart
│   └── docs_utils.jl       # Shared Julia utilities for docs
├── data/                   # Test input data
├── utils/                  # Utility scripts (test data download, SBOM, manifest)
└── sbom/                   # Software Bill of Materials
```

## Key Technologies & Dependencies

### Julia Stack (Core)
- **Julia ≥ 1.10**
- **NCDatasets** — Reading/writing NetCDF input and output files
- **BasicModelInterface** — BMI standard interface for model coupling
- **Graphs** — River network and graph-based routing
- **Polyester** — Lightweight multi-threading for computational loops
- **StaticArrays** — Stack-allocated small arrays for performance
- **CFTime** — Climate & Forecast time convention handling
- **TOML** — Configuration file parsing
- **Accessors / CompositionsBase** — Functional lenses for nested struct updates
- **Parameters** — Keyword-based struct constructors
- **FillArrays** — Memory-efficient filled arrays
- **Glob** — File path pattern matching
- **LoggingExtras / ProgressLogging / TerminalLoggers** — Structured logging and progress bars
- **OrderedCollections** — Insertion-ordered dictionaries
- **EnumX** — Enhanced enumerations
- **DelimitedFiles** — Reading delimited text data
- **PropertyDicts** — Dictionary-as-property access
- **Statistics** — Standard library statistics functions
- **Dates** — Date and time handling

### Custom internal systems (not external dependencies)

The codebase defines several internal utility systems instead of using external packages. Be aware of these before introducing new dependencies or reimplementing functionality:

- **Unit system** (`Wflow/src/units.jl`): A custom dimensional analysis system based on a `Unit` struct that represents physical units as products of powers (e.g. `Unit(; m = 1, s = -1)` for m/s). Supports SI conversion (`to_SI`, `from_SI`), arithmetic (`*`, `^`), string formatting (including BMI-standard output), and timestep-dependent units (`dt`). New units are added by extending the `Unit` struct and `to_SI_data`. Constants like `EMPTY_UNIT` and `ABSOLUTE_DEGREES` are predefined.

- **Standard name / parameter metadata system** (`Wflow/src/standard_name/`): An `OrderedDict`-based registry mapping string standard names to `ParameterMetadata` structs. Each `ParameterMetadata` records a lens (path into the model struct), unit, default value, fill value, type, description, and flags (e.g. `allow_dynamic_input`). Separate maps exist for SBM land hydrology, sediment, domain, and routing. Used throughout for netCDF I/O, TOML configuration binding, and documentation generation.

- **TOML configuration system** (`Wflow/src/config_structure.jl`, `config_init.jl`, `config_utils.jl`): The `Config` type wraps a parsed TOML as a hierarchy of `@kwdef mutable struct`s extending `AbstractConfigSection` (e.g. `TimeSection`, `ModelSection`, `InputSection`). TOML values are recursively converted to strongly-typed section structs via `init_config_section`, with automatic enum conversion via `EnumX`. The `InputEntry` struct handles three input modes: netCDF variable reference (with scale/offset/layer), uniform value, or external name.

- **NetCDF I/O layer** (`Wflow/src/io.jl`): `NCReader` manages input datasets and `Writer` manages output. `read_standardized` normalizes dimension ordering and direction (always increasing x/y). `ncread` combines config lookup, netCDF reading, metadata-driven defaults/fill values, and unit-aware type conversion into one call. The `param` helper navigates nested model structs by dot-separated path strings.

- **Graph and network utilities** (`Wflow/src/network.jl`, `routing/utils.jl`, `subdomains.jl`): `flowgraph` converts PCRaster-style LDD grids into `Graphs.DiGraph`. `NetworkLand` and `NetworkRiver` store index mappings between 1D internal and 2D external domains, staggered-grid edge connectivity (`EdgeConnectivity`), and subdomain decomposition for parallel kinematic-wave routing. `stream_order` computes Strahler order. `active_indices` derives forward/reverse index maps from subcatchment grids.

- **Threading utilities** (`Wflow/src/utils.jl`): `threaded_foreach` provides a unified parallel iteration abstraction—uses `Threads.@spawn` with chunked partitions when ≤ 8 threads, or `Polyester.@batch` for more. Used extensively in soil and routing update loops.

- **Numeric / hydrological helpers** (`Wflow/src/utils.jl`): `scurve` (sigmoid), `pow` (faster exponentiation via `exp(y*log(x))`), `tosecond` (Period → Float64 seconds), `bounded_divide`, `sin_slope`, `julian_day` (leap-day corrected), `lattometres` (lat/lon → meters), `set_layerthickness`, `svectorscopy` (Matrix → Vector{SVector}), `equal_size_vectors` (struct-of-arrays validation), and `water_table_change` (dynamic specific yield).

## Model architecture & simulation loop

The central type is `Model{R, L, M, T}`, parameterized on routing (`R <: Routing`), land model (`L <: AbstractLandModel`), mass balance (`M <: AbstractMassBalance`), and model type tag (`T <: AbstractModelType`). Fields: `config`, `domain`, `routing`, `land`, `mass_balance`, `clock`, `reader`, `writer`, `type`.

**Land vs Routing separation:**
- **Land models** (`AbstractLandModel`) compute vertical fluxes independently per cell. For hydrology, `LandHydrologySBM` contains: `atmospheric_forcing`, `vegetation_parameters`, `interception`, `snow`, `glacier`, `runoff`, `soil`, `demand`, `allocation`. For sediment, `SoilLossModel` is the land model.
- **Routing** (`Routing{O,R,S}`) bundles three horizontal flow components: `overland_flow`, `river_flow`, and `subsurface_flow`. Each can be a concrete type or a null `No*` type (e.g. `NoOverlandFlow`), avoiding conditional logic.

**Dispatch on model type tags:** Three singleton structs (`SbmModel`, `SbmGwfModel`, `SedimentModel`) are used as type parameters for dispatch. Constructors like `Model(config, ::SbmModel)` follow a consistent sequence: open static dataset → create `NCReader`/`Clock` → initialize `Domain` → initialize land model → initialize `Routing` → initialize mass balance → create `Writer` → close dataset → `set_states!` → return.

**Simulation loop** (`run_timestep!`):
1. `advance!(clock)` — increment time
2. `load_dynamic_input!(model)` — read forcing + cyclic parameters from NetCDF
3. `storage_prev!(model, mass_balance)` — snapshot storage for mass balance
4. `update_model!(model)` — dispatched on model type; updates land hydrology → exchanges fluxes to subsurface → updates subsurface flow → updates soil water → `surface_routing!` → `update_total_water_storage!`
5. `compute_mass_balance!(model, mass_balance)` — check conservation
6. `write_output(model)` — write NetCDF/CSV

**Immutable struct updates:** `@reset` from Accessors.jl is used extensively to "update" fields of immutable structs (creating copies with modified fields), especially during `Domain` initialization.

## Model types and routing

Three model types, set via `config.model.type`:
- **`sbm`**: Full hydrological model — SBM soil model + kinematic wave lateral subsurface flow + surface routing. Supports snow, glaciers, reservoirs, water demand/allocation.
- **`sbm_gwf`**: Like `sbm` but replaces lateral subsurface flow with a 2D unconfined aquifer groundwater flow model (`GroundwaterFlowModel`). Adds drain and constant-head boundary conditions.
- **`sediment`**: Soil erosion and sediment transport. Uses `SoilLossModel` as land model. No snow/glacier/subsurface. Mass balance is `NoMassBalance`.

Two routing types (`config.model.land_routing`, `config.model.river_routing`):
- **`kinematic_wave`** (default): Uses subdomains for parallel execution.
- **`local_inertial`**: Shallow water equations on an edge-based staggered grid. When both land and river use local inertial, a combined 2D routing method is dispatched.

Key configuration enums: `ModelType`, `RoutingType`, `VerticalConductivityProfile`, `GwfConductivityProfileType`, `RainfallErosionType`, `OverlandFlowErosionType`, `LandTransportType`, `RiverTransportType`.

## Domain and index conventions

`Domain` holds four sub-domains: `land` (`DomainLand`), `river` (`DomainRiver`), `reservoir` (`DomainReservoir`), `drain` (`DomainDrain`). Each has a Network struct and a Parameters struct.

**1D/2D index mapping:** The model internally works in 1D arrays over active cells only. `indices::Vector{CartesianIndex{2}}` maps from 1D internal domain to 2D external grid; `reverse_indices::Matrix{Int}` maps back (0 = inactive cell). No computation on inactive/masked cells.

**Network structs:**
- `NetworkLand`: `modelsize` (2D dims), `indices`, `reverse_indices`, `graph` (DiGraph), `order` (topological sort), `streamorder`, `upstream_nodes`, `order_of_subdomains`/`order_subdomain` (parallel kinematic wave), `river_indices` (land→river), `edge_indices` (staggered grid for local inertial).
- `NetworkRiver`: Similar plus `land_indices` (river→land), `reservoir_indices`, `nodes_at_edge`/`edges_at_node` (for local inertial), `pit_indices`.
- `NetworkReservoir`: `indices_outlet`, `indices_coverage`.

**Parameter structs:**
- `LandParameters`: `x_length`, `y_length`, `area`, `flow_width`, `surface_flow_width`, `flow_length`, `slope`, `flow_fraction_to_river`, `river_location`, `river_fraction`, `water_fraction`, `reservoir_outlet`/`coverage`.
- `RiverParameters`: `flow_width`, `flow_length`, `slope`, `cell_area`, `reservoir_outlet`/`coverage`.

## BMI interface

`Wflow/src/bmi.jl` implements the `BasicModelInterface.jl` protocol:
- `BMI.initialize(::Type{<:Model}, config_file)` → parses config, creates model, loads fixed forcing.
- `BMI.update(model)` → calls `run_timestep!`.
- `BMI.update_until(model, time)` → runs multiple timesteps to reach target time.
- `BMI.finalize(model)` → writes final state, resets clock, closes files.
- **Variable access:** `BMI.get_value_ptr` returns views into model arrays using the standard name metadata lens system. Variable names are CSDMS-style (e.g. `"river_water__volume_flow_rate"`).
- **Grid system:** 6 grid IDs: 0=reservoir, 1=drain, 2=river, 3–5=land. Variable-to-grid assignment is by string matching on the standard name.
- **`API` section in TOML config:** Lists which variables are exposed for BMI exchange.
- **Layered soil variables:** Names like `"soil_layer_2_water_unsaturated_zone__depth"` are parsed to extract the layer index and map to SVector storage.

## State management

- **Cold vs warm start:** Controlled by `config.model.cold_start__flag`. Cold start uses default values; warm start reads from `config.state.path_input`.
- **State identification:** States are identified by standard names and associated with tags (e.g. `:snow_state`, `:soil_state`, `:glacier_state`, `:vegetation_state`, `:reservoir_state`).
- **`extract_required_states(config)`:** Determines which states are needed based on model type and enabled features (snow, glacier, reservoir, floodplain, paddy).
- **`check_states(config)`:** Validates that the `[state.variables]` TOML section provides all required states. Throws `ArgumentError` if any are missing.
- **`set_states!(instate_path, model)`:** Reads from NetCDF, selects active indices, handles 3D (x,y,time) and 4D (x,y,layer,time) variables, converts to SVector for layered soil. Uses the metadata lens to set values.
- **Model-specific post-load:** `set_states!(model::AbstractModel{<:SbmModel})` computes derived variables (e.g. river storage from water depth, reservoir storage from water level).
- **State output:** Written by the `Writer` via `write_netcdf_timestep` at end of simulation or via `BMI.finalize`/`save_state`.

## Forcing and dynamic input

- **`AtmosphericForcing`** struct: `precipitation`, `potential_evaporation`, `temperature` (per cell).
- **Right-labeling convention:** Daily precipitation at `2000-01-02T00:00:00` is the total between 01-01 and 01-02.
- **`load_fixed_forcing!(model)`:** Sets time-invariant forcing values (uniform scalars from config) once at initialization.
- **`update_forcing!(model)`:** Reads time-varying NetCDF data at each timestep using backward-fill time interpolation.
- **`update_cyclic!(model)`:** Handles annually-repeating parameters (e.g. LAI). Reads from the static dataset at specific month-day combinations. Updates on first timestep and whenever the month-day matches an available time.
- **`load_dynamic_input!`:** Combines `update_forcing!` + `update_cyclic!`.
- **Reservoir mover pattern:** Precipitation and evaporation are averaged over reservoir coverage cells, set to zero in the land model, and stored separately in reservoir boundary conditions.
- **NCReader:** Supports multi-file datasets (`MFDataset` via glob patterns) and affine transforms (scale/offset) on input variables.

## Testing conventions

- **Framework:** Uses `TestItemRunner.jl` with `@testitem` macros (not standard `@testset`). Each `@testitem` is self-contained and can run independently. `runtests.jl` calls `@run_package_tests`.
- **Naming:** Unit tests are prefixed with `"unit: "` (e.g. `@testitem "unit: tosecond"`). Integration tests use descriptive names (e.g. `"Run SBM"`, `"BMI functions"`).
- **Test data:** Downloaded via `utils/download_test_data.jl` from GitHub releases (`wflow-artifacts`). Primarily Moselle and Piave catchment datasets in `Wflow/test/data/input/`.
- **Integration test pattern:** Initialize model → run 1–2 timesteps → assert specific values (using `≈` for floats) on model variables and output files.
- **State restart test:** Validates that a continuous run produces identical results to a cold-start + warm-start split (testing state save/load fidelity).
- **Code quality:** `Aqua.jl` tests check for ambiguities, unbound args, and other common issues.

## Julia code style
- Follow Julia community conventions
- Use multiple dispatch extensively
- Prefer immutable structs where possible
- Use `@kwdef` for struct definitions with defaults
- All code and inline documentation must be clear from a hydrological, mathematical and computer science perspective
- All new computational code must contain a docstring which specifies the units of all variables involved in the computation
- Single symbol iteration variables are not allowed. Use an informative name which is generally of the form `*_idx`, and see whether you can re-use such a variable name that is already used in similar code
- Keep the number of arguments to functions as low as possible. If a set of arguments is contained in a struct, only unpack the struct at the lowest level function where the individual arguments are actually needed. Do not create `NamedTuple`s to reduce the number of arguments
- Make the the argument types in signatures as precise as possible

## Julia performance considerations
- Avoid allocations during simulation
- Aim towards making it statically compilable
- Profile with `@profile` and `@benchmark`
- PrecompileTools.jl for reducing startup time
- Avoid type instabilities (use `@code_warntype`)

## Quarto documentation
- This repository uses Quarto to generate documentation
- This documentantion is mainly aimed at hydrologists
- The documentation contains Julia code blocks which are executed during the rendering of the documentation. These code blocks use the Wflow module and functionality defined at `docs/docs_utils.jl`. Please make sure that changes in these functionalities are handled properly in the Julia code blocks in the documentation
- For all Julia code in the documentation that is not executed but is there for clarification, if it is copied from an actual Julia script, make sure it stays up to date with that script

## Installation
- Installation is handled by pixi, see `pixi.toml`