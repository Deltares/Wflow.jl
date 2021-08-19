# Quick start

## Installation

Wflow is a [Julia](https://julialang.org/) package, that can be installed using Julia's
built-in package manager, Pkg. The rest of this section provides more detailed instructions,
such that new Julia users can also install Wflow. If you already know this, you can go ahead
to the next section.

First download and install the [current stable release of
Julia](https://julialang.org/downloads/#current_stable_release). Please see [platform
specific instructions](https://julialang.org/downloads/platform/) for further installation
instructions and if you have trouble installing Julia.

If you are new to Julia, it may be a good idea to check out the [Getting Started section of
the Julia Manual](https://docs.julialang.org/en/v1/manual/getting-started/). Links to other
learning resources can be found at
[julialang.org/learning](https://julialang.org/learning/).

Start an interactive session (also known as a read-eval-print loop or "REPL") by
double-clicking the Julia executable or running `julia` from the command line, in case you
have chosen to add it to your path. The default prompt will appear, in which Julia code can
be entered and run.

```julia-repl
julia> 1 + 2
3
```

Enter the Pkg REPL by pressing `]` from the Julia REPL. To get back to the Julia REPL, press
backspace or ^C. Just like in the Julia REPL, you can type `?` to get a help message, this
will list the available commands. To get more information on a specific command type `?`
followed by the command.

```julia-repl
pkg> ?

pkg> ?add
```

If you have not used Julia in a while, it can be a good idea to run `up` to update your
packages.
```julia-repl
pkg> up
```

If you go back to the Julia REPL, note that the same works, and is a convenient way to
consult documentation for specific functions or types, or other objects.

```
help?> cos
```

Now go back to the Pkg REPL and install Wflow using:

```julia-repl
pkg> add Wflow
```

This can take a while, especially the first time, since compatible dependencies are also
automatically looked up and installed from the Pkg General registry. It is also possible to
install Wflow from the `master` branch as follows:

```julia-repl
pkg> add https://github.com/Deltares/Wflow.jl
```

This command will track the `master` branch, and will update to the latest commit on that
branch when you run `update`, or simply `up`, in the Pkg REPL. The use of `add` will install
Wflow in you home directory under `.julia/packages/Wflow`. Note that packages installed
under `packages` by `add` are supposed to never be altered in that location, for Pkg and
it's automatic dependency handling to work well. If you want to make any changes to any of
the files in the Wflow directory, you need to do a development install. This can be done
using:

```julia-repl
pkg> dev https://github.com/Deltares/Wflow.jl
```

This will clone the git repository, put it under your home directory in `.julia/dev/Wflow`,
and add the Wflow package to your project environment. Note that to receive updates, you
have to pull in the latest changes yourself using `git pull`.

Finally, go back the Julia REPL and try to load Wflow:

```julia-repl
julia> using Wflow
```

The first time this will take longer as any package that is new or changed needs to be
precompiled first, to allow faster loading on subsequent uses.

You have now successfully installed Julia and Wflow. Before ending this section, we still
want to recommend a few tools that can make using and developing Julia code easier.

The first is `Revise.jl`. This package allows you to modify code and use the changes without
restarting Julia. Install it with `add Revise` from the Pkg REPL. Then create a file called
`.julia/config/startup.jl`, and put `using Revise` there. This will load Revise every time
you start a Julia session.

There is a section on editors and IDEs for Julia on https://julialang.org/, scroll down to
see it. We use and recommend Microsoft's free and open source [Visual Studio
Code](https://code.visualstudio.com/). When combined with the [Julia
extension](https://www.julia-vscode.org/) it provides a powerful and interactive development
experience.

## Usage

To run a simulation, first the required static and dynamic input data should be prepared in
NetCDF files. Next, we use a [TOML](https://github.com/toml-lang/toml) configuration file,
to fully specify a model simulation. This includes the model type, the start- and endtime,
the path to the forcing NetCDF files, as well as the mapping of parameters in the model to
names in the NetCDF file. See the section [Config and TOML](@ref) for more details.

As an example, we will run the Moselle example SBM model, of which the configuration file
can be seen here:
[sbm_simple.toml](https://github.com/Deltares/Wflow.jl/blob/master/test/sbm_simple.toml).
First we share some code that will download the required files to your current working
directory:

```julia
# urls to TOML and NetCDF of the Moselle example model
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.1/staticmaps.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/forcing-2000.nc"
toml_url = "https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_simple.toml"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "sbm_simple.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-moselle.nc"))
download(forcing, joinpath(datadir, "forcing-moselle.nc"))
download(toml_url, toml_path)
```

Now that we have our files in place, running a simulation is as simple as calling this
function:

```@docs
Wflow.run
```

As you can see, there are three different methods for this function. We will first make use
of the first one, where a path to a TOML file is passed as a `String`.

```julia
using Wflow
Wflow.run(toml_path)
```

This will parse the TOML file to create a `Wflow.Config`, use that to initialize a
`Wflow.Model`, and finally run the `Wflow.Model` for the specified duration.

```@docs
Wflow.Config
Wflow.Model
```

If you want to try out several different settings, without having to modify the TOML file
every time, you can create a `Wflow.Config` first, modify some settings, and then start the
simulation:

```julia
using Dates
config = Wflow.Config(toml_path)
config.endtime = DateTime("2000-01-03T00:00:00")
Wflow.run(config)
```

For even more control, you can initialize the model object yourself, and modify it directly,
or run a custom simulation loop. See the [`Wflow.run`](@ref) source for some
inspiration.

## [Multi-Threading](@id quickstart_multi_threading)

Wflow supports multi-threading execution of the wflow\_sbm model [SBM + Kinematic
wave](@ref) that uses the kinematic wave approach for river, overland and lateral subsurface
flow. Both the [SBM vertical concept](@ref) and the kinematic wave components of this model
can run on multiple threads. This functionality may also be useful for models that make
(partly) use of the kinematic wave as the [HBV model](@ref) and the wflow\_sbm model [SBM +
Groundwater flow](@ref).

For information on how to start Julia with multiple threads we refer to [How to start Julia
with multiple
threads](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads).
