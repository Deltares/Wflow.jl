# Build overview

## Directory tree
- `create_binaries`: folder containing code to compile the Wflow code
  - `create_binaries/create_app.jl`: compile and bundle wflow_cli with JuliaC
  - `create_binaries/wflow_cli.jl`: native executable entrypoint
  - `create_binaries/Project.toml`: project used to build wflow_cli
- `wflow_cli`: folder containing the code for the wflow_cli

The JuliaC build enables the workload in `Wflow/src/precompile.jl`, which runs representative
models during package precompilation to reduce wflow_cli startup latency.
