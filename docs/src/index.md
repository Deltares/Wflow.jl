```@meta
CurrentModule = Wflow
```

# About wflow

Wflow is Deltares’ solution for modelling hydrological processes, allowing users to account
for precipitation, interception, snow accumulation and melt, evapotranspiration, soil water,
surface water and groundwater recharge in a fully distributed environment. Successfully
applied worldwide for analyzing flood hazards, drought, climate change impacts and land use
changes, wflow is growing to be a leader in hydrology solutions. Wflow is conceived as a
framework, within which multiple distributed model concepts are available, which maximizes
the use of open earth observation data, making it the hydrological model of choice for data
scarce environments. Based on gridded topography, soil, land use and climate data, wflow
calculates all hydrological fluxes at any given grid cell in the model at a given time step.

Wflow was born out of the creation of Deltares in 2008, when a strategic review identified
the need for a distributed hydrological model to allow the simulation of flows at the
catchment scale. With the intention being to encourage greater scientific collaboration.
For this reason:

   * Wflow is free and open source software.
   * Wflow is easily coupled with other models and software applications.
   * Contribution to the wflow code development is encouraged.

From 2021 the [wflow code](https://github.com/Deltares/Wflow.jl) is distributed under the
[MIT License](https://github.com/Deltares/Wflow.jl/blob/master/LICENSE). Wflow is also
available as a [compiled executable](https://download.deltares.nl/en/download/wflow/) under
the Deltares terms and conditions. The wflow computational engine is built in the
[Julia](https://julialang.org/) language, a high-performance computing language and does not
have a graphical user interface, being designed for maximum user flexibility. Prior to 2021,
wflow was developed in PCRaster Python language and is [still
available](https://github.com/openstreams/wflow) but not actively developed.
