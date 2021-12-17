# [Tutorial models](@id sample_data)

## Overview of different model settings files

|  model  | TOML configuration file |
|:--------------- | ------------------|
| SBM + Kinematic wave | [sbm_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_config.toml) |      
| SBM + Groundwater flow | [sbm\_gwf\_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sbm_gwf_config.toml) | 
| HBV model | [hbv_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/hbv_config.toml) | 
| wflow_sediment | [sediment_config.toml](https://raw.githubusercontent.com/Deltares/Wflow.jl/master/test/sediment_config.toml) |



## [wflow\_sbm](@id wflow_sbm_data)

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

## wflow\_sbm + Groundwater flow
```julia
# urls to TOML and NetCDF of the Moselle example model
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
# urls to TOML and NetCDF of the Moselle example model
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

## wflow\_sediment
```julia
# urls to TOML and NetCDF of the Moselle example model
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