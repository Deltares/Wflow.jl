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
│   starttime = 2000-01-01T00:00:00
│   Δt = 86400 seconds
└   endtime = 2000-02-01T00:00:00

Progress: 100%|██████████████████████████████████████████████████| Time: 0:00:27
```
