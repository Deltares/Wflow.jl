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
