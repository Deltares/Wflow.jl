# Trimmed Wflow binary (proof of concept)

Build: `pixi run -e trim build-trim-poc`.
Run: `build/trim_poc/out/wflow_trim/bin/wflow_trim path/to/config.toml`.

The supported model options are fixed at compile time by a build profile
(`PROFILE` in `wflow_trim/src/wflow_trim.jl`, overridable with `WFLOW_TRIM_PROFILE`
at build time). Verified to give results identical to `release/v1.1`:

- default profile (SBM, kinematic wave, no snow/reservoirs): `Wflow/test/sbm_trim_poc.toml`
- `(snow = true, reservoir = true)`: `sbm_config.toml`, cold and warm start
- plus local inertial river and land routing: `sbm_river-land-local-inertial_config.toml`

Not yet trim compatible: water demand, sbm_gwf, sediment, non-standard calendars, BMI.

## Measurements

Julia 1.13.1, JuliaC 0.3.10, 32 threads, Linux, based on `release/v1.1` before #1118.

| | untrimmed `wflow_cli` (release/v1.1) | trimmed `wflow_trim` |
|---|---|---|
| build time | 313 s | 76 s (40 s with native CPU target) |
| build peak memory | 17.5 GB | 3.5 GB |
| executable | 724 MB | 37 MB (19 MB native) |
| bundle (zstd -10) | 1.1 GB (245 MB) | 215 MB (56 MB) |
| 9 day POC run, wall time / memory | 3.1 s / 760 MB | 1.1 s / 210 MB |
| one year POC run, wall time | 63 s | 56 s |

Both builds use multi CPU targets. The untrimmed build supports all model types and runs a
precompile workload (needs the test data); the trimmed one supports one build profile.

In normal Julia (`using Wflow`, same configs) this branch vs release/v1.1: first 9 day run
including compilation 57 vs 86 s, one year run 54 vs 61 s, allocations 4.7 vs 7.0 GB. So
the runtime gain comes from the type stability changes, not from trimming itself.

## Logging (`wflow_trim/src/trim_logging.jl`)

JuliaC disables the logging macros under trim by overriding
`Base.CoreLogging.current_logger_for_env(level, group, _module) = nothing`, and the macros
call the logger via `@nospecialize` and `@invokelatest`. The app adds more specific
methods for one concrete `TrimLogger` (`current_logger_for_env(::LogLevel, ::Symbol,
::Module)`, `handle_message_nothrow`, `logging_error`), which makes `@info` and friends
dispatch statically. It writes to the terminal and to `path_log` in the `ConsoleLogger`
layout, and respects `loglevel` and `silent`. This is type piracy that is only acceptable
in the binary entry point.

Enabling logging makes every log statement reachable, so it must be trim compatible too:
no `show` of arrays or Dates periods in messages (`vector_string`, `format_duration`), and
concretely inferred keyword argument values. Log statements of dependencies stay disabled
(`is_logged_module`).

## Working on trim compatibility

Workflow:

- Build with another profile, e.g.
  `WFLOW_TRIM_PROFILE='(snow = true, reservoir = true)' pixi run -e trim build-trim-poc`.
- On failure, the build prints `Trim verify finished with N errors` and a stack trace per
  `Verifier error`. The first frame in Wflow shows the culprit; many errors usually share a
  few causes, so fix those and rebuild (about 1 minute).
- To check results, run the same TOML with normal Julia (`Wflow.run`) and compare the output
  files, e.g. with NCDatasets.
- To check that the bundle is relocatable, run it with an empty `HOME`; setting
  `JULIA_DEPOT_PATH` gives false failures.

Patterns that work under `--trim=safe`:

- Values from a parsed TOML are `Any`. Dispatch on the target type only, mark the value
  `@nospecialize` and branch with `isa` (see `convert_value` in `config_init.jl`). A call
  with an `Any` argument only resolves if the callee's argument is `@nospecialize`.
- Config values that select types go through `option(type, Val(:name), value)`, so that the
  build profile turns them into constants.
- Struct fields must be concrete: add type parameters, or use small `Union`s.
- Add typeasserts where inference loses the type. The metadata lookup in `ncread` uses
  `Base.@constprop :aggressive` and `Base.@assume_effects :foldable`, but it doesn't always
  fold.
- Keyword arguments, including log keyword arguments, need concretely inferred values. Use
  `something(x)` or a local variable after the `isnothing` check.
- Use a `Task[]` vector instead of `@sync`, which waits on a `Channel{Any}`.
- `@generated` functions must be defined after the globals they read
  (`standard_name_switch.jl` is included after the standard name maps).
- Dynamic model paths (`param(model, "a.b.c")`) only work outside a build profile, through
  `with_model_field`.

Known not to trim: NCDatasets and CFTime (time decoding and formatting), `string(::DateTime)`,
`Dates.canonicalize`, `show` of arrays, `TOML.print`, DelimitedFiles, Glob, `readlines`,
TerminalLoggers and LoggingExtras, and the error path of `Statistics.median!`.

Known test failures in normal Julia:

- The unit tests in `Wflow/test/io.jl` that use the old NCDatasets and CFTime internals
  (`get_at`, `cftime`, `timecycles`, `read_x_axis`, ...) are not ported yet.
- Tests that need fresh data (`dvc pull`) fail with stale local test data.

## Findings

- Small `Union` fields pass the verifier. Abstractly typed fields only work for calls
  through world splitting (at most `max_methods` matching methods, default 3), not for
  field access through them.
- Julia 1.14 nightly: the POC builds with one workaround and gives the same output. Fewer
  libraries are loaded at startup, but the bundle is still ~200 MB, as JuliaC cannot link
  the new static `libjulia-internal.a` yet.
