# Trimmed Wflow binary (proof of concept)

Build: `pixi run -e trim build-trim-poc` (or `build/trim_poc/run_build.sh wflow_trim <log>`).
Output: `build/trim_poc/out/wflow_trim/bin/wflow_trim path/to/config.toml`.

The supported model options are fixed at compile time by a build profile
(`PROFILE` in `wflow_trim/src/wflow_trim.jl`, overridable with `WFLOW_TRIM_PROFILE`
at build time). Verified to give results identical to `release/v1.1`:

- default profile (SBM, kinematic wave, no snow/reservoirs): `Wflow/test/sbm_trim_poc.toml`
- `(snow = true, reservoir = true)`: `sbm_config.toml`, cold and warm start
- plus local inertial river and land routing: `sbm_river-land-local-inertial_config.toml`

Not yet trim compatible: water demand, sbm_gwf, sediment, non-standard calendars, BMI.
`baseline/`, `full/` and `sbm_all/` measure remaining verifier errors; `probe_*`
test dependencies in isolation.

## Logging (`wflow_trim/src/trim_logging.jl`, `probe_log/`)

JuliaC disables the logging macros under trim by overriding
`Base.CoreLogging.current_logger_for_env(level, group, _module) = nothing`. The macros
themselves call the logger via `@nospecialize` and `@invokelatest`, so a custom logger
type alone is not enough. The app adds more specific methods for one concrete
`TrimLogger` (`current_logger_for_env(::LogLevel, ::Symbol, ::Module)`,
`handle_message_nothrow`, `logging_error`), which makes `@info` and friends dispatch
statically. It writes to the terminal and to `path_log`, respects `loglevel` and
`silent`, and uses the `ConsoleLogger` layout; the log file of the POC run matches
`release/v1.1` apart from `DateTime` instead of `DateTimeStandard` values and the
2D output index logged as `x`, `y`. Works on 1.13.1 and 1.14 nightly; it
is type piracy that is only acceptable in the binary entry point.

Enabling logging makes the message and keyword argument code of every log statement
reachable, so it must be trim compatible too:

- No `show` of arrays or Dates periods in messages (`vector_string`, `format_duration`).
- Keyword argument values must be concretely inferred, a `Union{Nothing, Int}` field
  gives an unresolved `kwcall` (use a local variable after the `isnothing` check).
- Keyword values are printed for common types only (strings, numbers, `DateTime`,
  exceptions), others by type name.
- Log statements of dependencies are kept disabled (`is_logged_module`), e.g. the
  `@debug` in `ThreadingUtilities.__init__` builds its message with an untrimmable
  `string` call.

Upstream: Julia PR #55166 ("make logging work with trimming", normalizing messages and
keyword arguments to strings) was closed unmerged in 2024; no follow-up is open.

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
- `with_logger` now passes the verifier, but JuliaC still disables logging under trim
  (worked around with `trim_logging.jl`, see above).
