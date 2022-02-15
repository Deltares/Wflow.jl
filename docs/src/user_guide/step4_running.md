# [Step 4: Running a simulation](@id run_wflow)

## Using Julia

Once you installed Julia and Wflow.jl, a simulation can be started from the command line
as follows:

```
julia -e 'using Wflow; Wflow.run()' path/to/config.toml
``` 

Furthermore, it is possible to write a Julia script to run a simulation. Example data to 
explore how this works can be found [here](@ref sample_data).

```julia
using Wflow
Wflow.run(toml_path)
```

Julia can also be used to modify settings after reading the settings file. In the example
below, we show how to adjust the end date of the simulation.

```julia 
using Dates
config = Wflow.Config(toml_path)
config.endtime = DateTime("2000-01-03T00:00:00")
Wflow.run(config)
```

## [Using the command line interface](@id cli)

If you don't need the extra features of using Wflow as a library, but just want to run
simulations, the command line interface makes it easier to do so. It consists of a single
executable, `wflow_cli` that accepts a single argument, the path to a TOML configuration
file.

Binaries of `wflow_cli` can be downloaded from our website
[download.deltares.nl](https://download.deltares.nl/en/download/wflow/), and are currently
available for Windows.

After installing you can see three folders in the installation directory. It is only the
`bin/wflow_cli` that is directly used. All three folders need to stay together however.
The share folder contains TOML files with more information about the build.

```
bin\wflow_cli
lib
share
```

Simply running `wflow_cli` with no arguments will give the following message:

```
Usage: wflow_cli 'path/to/config.toml'
```

When starting a run, you will see basic run information on the screen, as well as a progress
bar, that gives an estimate of how much time is needed to finish the simulation:

```
┌ Info: Run information
│   model_type = "sbm"
│   starttime = CFTime.DateTimeStandard(2000-01-01T00:00:00)
│   Δt = 86400 seconds
│   endtime = CFTime.DateTimeStandard(2000-12-31T00:00:00)
└   nthreads() = 4

Progress: 100%|██████████████████████████████████████████████████| Time: 0:00:27
```
