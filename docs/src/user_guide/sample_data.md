# [Example models](@id sample_data)

For each Wflow Model a test model is available that can help to understand the data
requirements and the usage of each Model. The TOML configuration file per available model
are listed in the Table below:

|  model  | TOML configuration file |
|:--------------- | ------------------|
| wflow\_sbm + kinematic wave | [sbm_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_config.toml) |
| wflow\_sbm + groundwater flow | [sbm\_gwf\_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_gwf_config.toml) |
| wflow\_hbv | [hbv_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/hbv_config.toml) |
| wflow\_flextopo | [flextopo_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/flextopo_config.toml) |
| wflow_sediment | [sediment_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sediment_config.toml) |

The associated Model files (input static, forcing and state files) can easily be downloaded
and for this we share the following Julia code (per Model) that downloads the required files
to your current working directory. For running these test model see also [Usage](@ref run_wflow)
and [Command Line Interface](@ref cli).

## [wflow\_sbm + kinematic wave](@id wflow_sbm_data)
```julia
# urls to TOML and netCDF of the Moselle example model
toml_url = "https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_config.toml"
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.7/staticmaps-moselle.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.6/forcing-moselle.nc"
instates = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.6/instates-moselle.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "sbm_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-moselle.nc"))
download(forcing, joinpath(datadir, "forcing-moselle.nc"))
download(instates, joinpath(datadir, "instates-moselle.nc"))
download(toml_url, toml_path)
```

## wflow\_sbm + groundwater flow
```julia
# urls to TOML and netCDF of the Moselle example model
toml_url = "https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_gwf_config.toml"
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.2/staticmaps-sbm-groundwater.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.1/forcing-sbm-groundwater.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "sbm_gwf_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-sbm-groundwater.nc"))
download(forcing, joinpath(datadir, "forcing-sbm-groundwater.nc"))
download(toml_url, toml_path)
```

## wflow\_hbv
```julia
# urls to TOML and netCDF of the Moselle example model
toml_url = "https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/hbv_config.toml"
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.1/staticmaps-lahn.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/forcing-lahn.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "hbv_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-lahn.nc"))
download(forcing, joinpath(datadir, "forcing-lahn.nc"))
download(toml_url, toml_path)
```

## [wflow\_flextopo](@id wflow_flextopo_data)
```julia
# urls to TOML and netCDF of the Meuse example model
toml_url = "https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/flextopo_config.toml"
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.8/staticmaps_flex_meuse.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.8/forcing_meuse.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "flextopo_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps_flex_meuse.nc"))
download(forcing, joinpath(datadir, "forcing_meuse.nc"))
download(toml_url, toml_path)
```

## wflow\_sediment
```julia
# urls to TOML and netCDF of the Moselle example model
toml_url = "https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sediment_config.toml"
staticmaps = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.3/staticmaps-moselle-sed.nc"
forcing = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.3/forcing-moselle-sed.nc"
instates = "https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/instates-moselle-sed.nc"

# create a "data" directory in the current directory
datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
toml_path = joinpath(@__DIR__, "sediment_config.toml")

# download resources to current and data dirs
download(staticmaps, joinpath(datadir, "staticmaps-moselle-sed.nc"))
download(forcing, joinpath(datadir, "forcing-moselle-sed.nc"))
download(instates, joinpath(datadir, "instates-moselle-sed.nc"))
download(toml_url, toml_path)
```
