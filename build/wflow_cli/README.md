# wflow_cli

This is a [Julia](https://julialang.org/) project that uses the
[Wflow.jl](https://github.com/Deltares/Wflow.jl) Julia package, puts a simple command line
interface (cli) on top, and packages this into a standalone application using
[PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl).

This enables using Wflow.jl without having to install Julia, and thus makes it more
convenient to use in certain settings where installation must be simple and no interactive
Julia session is needed.

If you have installed Julia and Wflow.jl, a simulation can also be started from the command
line as follows:

```
julia -e 'using Wflow; Wflow.run()' path/to/config.toml
```

With a wflow_cli build this becomes:

```
wflow_cli path/to/config.toml
```

## Directory tree
- `Setup/`: used to create MSI installers for Windows
- `src/`: source code of the wflow_cli app, that wraps `Wflow.run` function
- `Project.toml`: Project for the wflow_cli app that needs to be built
