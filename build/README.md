# Build overview

## Directory tree
- `create_binaries`: folder containing code to compile the Wflow code
  - `create_binaries/create_app.jl`: compiles and bundles wflow_cli with JuliaC
- `wflow_cli`: project containing the wflow_cli source and native executable entrypoint

The wflow_cli package runs representative models during package precompilation to reduce startup
latency. The workload is build-specific because it uses downloaded test data; installing or
developing Wflow does not run it.
