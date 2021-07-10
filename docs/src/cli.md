# Command Line Interface

If you don't need the extra features of using Wflow as a library, but just want to run
simulations, the command line interface makes it easier to do so. It consists of a single
executable, `wflow_cli` that accepts a single argument, the path to a TOML configuration
file.

Binaries of `wflow_cli` can be downloaded from our website
[download.deltares.nl](https://download.deltares.nl/en/download/wflow/), and are currently
available for Windows and Linux.

After installing you can see two folders in the installation directory. It is only the
`bin/wflow_cli` that is used. The artifacts folder contains binary dependencies such as
NetCDF.

```
artifacts\
bin\wflow_cli
```

Simply running `wflow_cli` with no arguments will give the following message:

```
Usage: wflow_cli 'path/to/config.toml'
```

When starting a run, you will see basic run information on the screen, as well as a progress
bar, that gives an estimate of how much time is needed to finish the simulaion:

```
┌ Info: Run information
│   model_type = "sbm"
│   starttime = CFTime.DateTimeStandard(2000-01-01T00:00:00)
│   Δt = 86400 seconds
│   endtime = CFTime.DateTimeStandard(2000-12-31T00:00:00)
└   nthreads() = 4

Progress: 100%|██████████████████████████████████████████████████| Time: 0:00:27
```

To try out a simple test model, you can download
[`sbm_moselle_config_data.zip`](https://github.com/visr/wflow-artifacts/releases/download/v0.2.1/sbm_moselle_config_data.zip),
which includes both a TOML configuration file as well as the NetCDF input data, for a one
year simulation of the Moselle catchment.

## [Multi-Threading](@id cli_multi_threading)

As explained in the [quick start section on multi-threading](@ref quickstart_multi_threading), we need to start julia
with multiple threads to make use of this speedup. For `wflow_cli`, the only way to do this
is by setting the `JULIA_NUM_THREADS` environment variable, as explained in [these julia
docs](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads).

When a model run starts, among the run information the number of threads that are used is
printed, so `nthreads() = 4` means 4 threads are used, because `JULIA_NUM_THREADS` has been
set to 4.
