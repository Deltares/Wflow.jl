# Build overview

## Directory tree
- `create_binaries`: folder containing code to compile the Wflow code
  - `create_binaries/create_app.jl`: compiles and bundles wflow_cli with JuliaC
- `wflow_cli`: project containing the wflow_cli source and native executable entrypoint

The wflow_cli project opts into the workload in `Wflow/src/precompile.jl`, which runs
representative models during package precompilation to reduce startup latency. Regular Wflow
development does not run this workload.
