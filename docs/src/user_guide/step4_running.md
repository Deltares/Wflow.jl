# [Step 4: Running a simulation](@id run_wflow)

## Using Julia

Below shows a brief example on how to run wflow using Julia. Example data to explore how this works can be found [here](@ref sample_data).

```julia
using Wflow
Wflow.run(toml_path)
```

Julia can also be used to modify settings after reading the settings file. In the example below, we show how to adjust the end date of the simulation.

```julia 
using Dates
config = Wflow.Config(toml_path)
config.endtime = DateTime("2000-01-03T00:00:00")
Wflow.run(config)
```

## [Using the command line interface](@id cli)

Simply running `wflow_cli` in the command line with no arguments will give the following 
message:

```
Usage: wflow_cli 'path/to/config.toml'
```

To try out a simple test model, you can download the following example data:
[`sbm_moselle_config_data.zip`](https://github.com/visr/wflow-artifacts/releases/download/v0.2.1/sbm_moselle_config_data.zip).
This contains the following files, with the input data in the `data` folder, and the model
settings in the `.toml` file:

```
data\
    forcing-moselle.nc
    staticmaps-moselle.nc
sbm_config.toml
```

The TOML configuration file is prepared with the settings required to perform a one
year simulation of the Moselle catchment. To run the model, execute the following command. 
This will give the following message with information on the run, while the model is 
simulating the catchment

```
wflow_cli "path/to/sbm_moselle/sbm_config.tml"

┌ Info: Run information
│   model_type = "sbm"
│   starttime = CFTime.DateTimeStandard(2000-01-01T00:00:00)
│   Δt = 86400 seconds
└   endtime = CFTime.DateTimeStandard(2000-12-31T00:00:00)

Progress: 100%|██████████████████████████████████████████████████| Time: 0:01:07
```
