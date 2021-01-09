```@meta
CurrentModule = Wflow
```

# Wflow

Wflow is a [Julia](https://julialang.org/) package that provides a hydrological modeling
framework, as well as several different vertical and lateral concepts that can be used to
run hydrological simulations. It is a continuation of the work available
[here](https://github.com/openstreams/wflow).

You can use Wflow in the following three ways.

1. From Julia as a library. The [Quick start](@ref) describes how to
   install it, prepare an example model, and run a simulation.
2. With the `wflow_cli` [Command Line Interface](@ref). This way there is no
   need to install Julia, you can download a binary and run.
3. Through the [Basic Model Interface](@ref), which is a standardized interface to run the
   model that also makes it easy to couple other BMI enabled models to it.

Wflow can also be linked to the Flood forecasting system
[Delft-FEWS](https://oss.deltares.nl/web/delft-fews/). This is described in [Run from
Delft-FEWS](@ref)

See the side bar on the left to navigate to other sections of the documentation.