# Trimmed Wflow binary (proof of concept)

Build: `pixi run -e trim build-trim-poc` (or `build/trim_poc/run_build.sh wflow_trim <log>`).
Output: `build/trim_poc/out/wflow_trim/bin/wflow_trim path/to/config.toml`.

The supported model options are fixed at compile time by a build profile
(`PROFILE` in `wflow_trim/src/wflow_trim.jl`, overridable with `WFLOW_TRIM_PROFILE`
at build time). Verified to give results identical to `release/v1.1`:

- default profile (SBM, kinematic wave, no snow/reservoirs): `Wflow/test/sbm_trim_poc.toml`
- `(snow = true, reservoir = true)`: `sbm_config.toml`, cold and warm start
- plus local inertial river and land routing: `sbm_river-land-local-inertial_config.toml`

Not yet trim compatible: water demand, sbm_gwf, sediment, non-standard calendars,
logging (JuliaC disables logging macros in trimmed builds), BMI.
`baseline/`, `full/` and `sbm_all/` measure remaining verifier errors; `probe_*`
test dependencies in isolation.

## Findings on dispatch (juliac 1.13.1, `probe_dispatch/`)

- Small `Union` fields pass the verifier (tested with 2 and 5 members).
- Abstractly typed fields pass only through world splitting: a call resolves if at most
  `max_methods` methods match (default 3, can be raised per function with
  `Base.Experimental.@max_methods n function f end`, tested up to 5). Field access through an
  abstract field (`x.snow.variables`) is not typed, so only calls benefit.
- A manual `isa` split at the call site always works.

## Julia 1.14 (nightly 1.14.0-DEV.3678, JuliaC 0.3.11)

- The POC builds and gives the same output, after one workaround (`print` of a
  `Union{Nothing, String}` was reported unresolved, see `config_utils.jl`).
- OpenBLAS, SuiteSparse and libgit2 are no longer loaded at startup (98 vs 153 MB of loaded
  libraries). The bundle is still ~200 MB; JuliaC 0.3.11 cannot link the new static
  `libjulia-internal.a` yet. The executable is larger (37 vs 17 MB).
- `with_logger` now passes the verifier, but JuliaC still disables logging under trim.
