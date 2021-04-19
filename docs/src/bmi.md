# Basic Model Interface

## Introduction
The [Community Surface Dynamics Modeling System](https://csdms.colorado.edu/wiki/Main_Page)
(CSMDS) has developed the Basic Model Interface (BMI). BMI consists of a set of standard
control and query functions that can be added by a developer to the model code and makes a
model both easier to learn and easier to couple with other software elements. 

For more information see also: <http://csdms.colorado.edu/wiki/BMI_Description>

CSDMS provides specifications for the languages C, C++, Fortran and Python. Wflow, written
in the [Julia programming language](https://julialang.org/), makes use of the following
[Julia specification](https://github.com/Deltares/BasicModelInterface.jl), based on BMI 2.0
version.

For the BMI implementation of Wflow all grids are defined as [unstructured
grids](https://bmi-spec.readthedocs.io/en/latest/model_grids.html#unstructured-grids). While
the input (forcing and model parameters) is structured (uniform rectilinear), internally
wflow works with one dimensional arrays based on the active grid cells of the 2D model
domain.

## Configuration
The variables that Wflow can exchange through BMI are based on the different model
components and these components should be listed under the `API` section of the TOML
configuration file of the model type. Below an example of this `API` section, that lists the
`vertical` component and different `lateral` components:

```toml
[API]
components = [
  "vertical",
  "lateral.subsurface",
  "lateral.land",
  "lateral.river",
  "lateral.river.reservoir"
]
```

See also:
```@docs
Wflow.BMI.initialize
Wflow.BMI.get_input_var_names
```

## Couple to a groundwater model
For the coupling of wflow\_sbm (SBM + kinematic wave) with a groundwater model (e.g.
MODFLOW) it is possible to run:
- wflow\_sbm in parts from the BMI, and 
- to switch off the lateral subsurface flow component of wflow\_sbm. 

The lateral subsurface component of wflow\_sbm is not initialized by Wflow when the
`[input.lateral.subsurface]` part of the TOML file is not included. Then from the BMI it is
possible to run first the recharge part of SBM:

```julia
model = BMI.update(model, run="sbm_until_recharge")
```
and to exchange recharge and for example river waterlevels to the groundwater model. After
the groundwater model update, and the exchange of groundwater head and for example drain and
river flux to wflow\_sbm, the SBM part that mainly determines exfiltration of water from the
unsaturated store, and the kinematic wave for river - and overland flow can be run as
follows: 

```julia
model = BMI.update(model, run="sbm_after_subsurfaceflow")
```

See also:
```@docs
Wflow.BMI.update
```