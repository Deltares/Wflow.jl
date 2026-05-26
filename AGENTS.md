# Copilot instructions

Wflow — Julia hydrological modeling framework (Julia ≥ 1.10). Installation via `pixi.toml`.

## Repository layout
- `Wflow/src/` — core package source
- `Wflow/test/` — tests (`TestItemRunner.jl`, `@testitem` macros)
- `build/` — binary compilation scripts
- `server/` — ZMQ-based BMI server
- `docs/` — Quarto documentation (targets hydrologists); Julia code blocks execute using `docs/docs_utils.jl`
- `data/`, `utils/` — test data and utility scripts

## Key dependencies
NCDatasets, BasicModelInterface, Graphs, Polyester, StaticArrays, CFTime, Accessors, Parameters, EnumX, OrderedCollections

## Internal systems (do not reimplement)
- **Unit system** (`units.jl`): `Unit` struct for dimensional analysis, `to_SI`/`from_SI`, timestep-dependent `dt` units
- **Standard name metadata** (`standard_name/`): `OrderedDict` mapping standard names → `ParameterMetadata` (lens, unit, default, type, flags). Used for NetCDF I/O, TOML binding, docs
- **Config system** (`config_structure.jl`, `config_init.jl`, `config_utils.jl`): `Config` wraps TOML as typed `AbstractConfigSection` structs. `InputEntry` handles netCDF ref / uniform value / external name
- **NetCDF I/O** (`io.jl`): `NCReader`/`Writer`; `ncread` combines config lookup + reading + defaults + unit conversion
- **Network/graph** (`network.jl`, `routing/utils.jl`, `subdomains.jl`): `flowgraph` converts LDD → DiGraph; `NetworkLand`/`NetworkRiver` store 1D↔2D index maps, edge connectivity, subdomain decomposition
- **Threading** (`utils.jl`): `threaded_foreach` — `Threads.@spawn` (≤8 threads) or `Polyester.@batch`
- **Numeric helpers** (`utils.jl`): `scurve`, `pow`, `tosecond`, `bounded_divide`, `lattometres`, `svectorscopy`, etc.

## Architecture
- Central type: `Model{R, L, M, T}` — routing, land model, mass balance, model type tag
- **Land models** (`AbstractLandModel`): vertical per-cell fluxes. `LandHydrologySBM` (hydrology), `SoilLossModel` (sediment)
- **Routing** (`Routing{O,R,S}`): `overland_flow`, `river_flow`, `subsurface_flow` — each concrete or null (`No*` type)
- **Model type tags** (dispatch singletons): `SbmModel`, `SbmGwfModel`, `SedimentModel`
- **Immutable struct updates**: use `@reset` from Accessors.jl

## Simulation loop (`run_timestep!`)
1. `advance!(clock)`
2. `load_dynamic_input!`
3. `storage_prev!`
4. `update_model!`
5. `compute_mass_balance!`
6. `write_output`

## Model types (`config.model.type`)
- `sbm`: SBM soil + kinematic wave subsurface + surface routing (snow, glaciers, reservoirs, demand)
- `sbm_gwf`: like `sbm` but 2D groundwater flow replaces lateral subsurface
- `sediment`: soil erosion/sediment transport, `NoMassBalance`

## Routing (`config.model.land_routing`, `config.model.river_routing`)
- `kinematic_wave` (default): parallel via subdomains
- `local_inertial`: shallow water equations, staggered grid; combined 2D when both land+river use it

## Domain & indexing
- `Domain` has `land`, `river`, `reservoir`, `drain` sub-domains (each: Network + Parameters)
- Internal 1D arrays over active cells only; `indices` (1D→2D), `reverse_indices` (2D→1D, 0=inactive)

## BMI (`bmi.jl`)
- `BMI.initialize`, `update`, `update_until`, `finalize`
- `get_value_ptr` uses standard name lens system; CSDMS-style names
- 6 grid IDs: 0=reservoir, 1=drain, 2=river, 3–5=land
- `API` TOML section lists exposed variables
- Layered soil: name pattern `"soil_layer_N_..."` → SVector index

## States
- Cold/warm start: `config.model.cold_start__flag`; warm reads from `config.state.path_input`
- `extract_required_states(config)`, `check_states(config)`, `set_states!(path, model)`
- 3D (x,y,time) and 4D (x,y,layer,time) NetCDF variables; layered → SVector

## Forcing
- `AtmosphericForcing`: precipitation, potential_evaporation, temperature
- Right-labeling: timestamp marks end of accumulation period
- `load_dynamic_input!` = `update_forcing!` + `update_cyclic!`
- Reservoir: precip/evap averaged over coverage cells, zeroed in land model

## Running Wflow

All commands use `pixi run` which activates the correct Julia environment.

**Loading/compiling the package** (verifies no errors):
```
pixi run julia --project=Wflow -e 'using Wflow; println(\"Wflow loaded successfully\")'
```

**Shell quoting (Windows PowerShell)**: Use outer single quotes with inner escaped
double quotes (`\"`). PowerShell mangles nested double quotes; triple-quotes and
`raw""` strings do not work reliably.

## Testing
- To run all tests:

```
pixi run julia --project=Wflow -e 'using Pkg; Pkg.test()'
```

- To run a subset of test items filtered by name (e.g. unit tests only):

```
pixi run julia --project=Wflow -e 'using Pkg; Pkg.test(test_args=[\"unit\"])'
```

- `TestItemRunner.jl` with `@testitem` (not `@testset`); unit tests prefixed `"unit: "`
- Test data: `utils/download_test_data.jl` → Moselle/Piave datasets
- Integration: init model → run 1–2 steps → assert values with `≈`
- State restart test: continuous run == cold-start + warm-start split
- `Aqua.jl` for code quality checks

## Code style (mandatory)
- Multiple dispatch; prefer immutable structs; `@kwdef` for defaults
- Precise argument types in signatures
- Docstrings on all computational code specifying units of all variables
- No single-symbol iteration variables — use descriptive `*_idx` names
- Minimize function arguments; pass structs, unpack only at lowest level; no `NamedTuple` wrappers
- Code must be clear from hydrological, mathematical, and CS perspectives

## Performance (mandatory)
- Zero allocations during simulation loops
- No type instabilities (`@code_warntype`)
- Target static compatibility
- PrecompileTools.jl for startup

## Documentation
- Quarto-based; Julia code blocks execute during render using `docs/docs_utils.jl`
- Keep non-executed example code in sync with actual source files
