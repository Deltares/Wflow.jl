# Build overview

## Directory tree
- `create_binaries`: folder containing code to compile the Wflow code
  - `create_binaries/create_app.jl`: run PackageCompiler
  - `create_binaries/precompile.jl`: run models to compile code, making wflow_cli startup faster
  - `create_binaries/Project.toml`: project used to build wflow_cli
- `wflow_cli`: folder containing the code for the wflow_cli
