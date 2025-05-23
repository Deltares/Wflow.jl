# Wflow

[![Build Status](https://github.com/Deltares/Wflow.jl/workflows/CI/badge.svg)](https://github.com/Deltares/Wflow.jl/actions)
[![Coverage](https://codecov.io/gh/Deltares/Wflow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Deltares/Wflow.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://deltares.github.io/Wflow.jl/dev)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://deltares.github.io/Wflow.jl/stable)
[![DOI](https://zenodo.org/badge/246787232.svg)](https://zenodo.org/badge/latestdoi/246787232)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Wflow is a [Julia](https://julialang.org/) package that provides a hydrological modeling
framework, as well as several different vertical and lateral concepts that can be used to
run hydrological simulations. It is a continuation of the work available
[here](https://github.com/openstreams/wflow).

## Documentation
See [stable](https://deltares.github.io/Wflow.jl/stable) or
[dev](https://deltares.github.io/Wflow.jl/dev) for the documentation.

## Installation
For the installation enter the Pkg REPL by pressing `]` from the Julia REPL, and then
```julia-repl
pkg> add Wflow
```
A more detailed description, including a development install, is available
[here](https://deltares.github.io/Wflow.jl/dev/getting_started/install).

## Contributions and reporting issues
We welcome reporting of issues [here](https://github.com/Deltares/Wflow.jl/issues). Please
provide a minimum working example so we are able to reproduce the issue. Furthermore, we
welcome contributions. We follow the [ColPrac guide for collaborative
practices](https://github.com/SciML/ColPrac). New contributors should make sure to read that
guide.

## Citing
For publications, please cite the following paper introducing Wflow.jl and describing the
wflow\_sbm concept, together with some case studies:

van Verseveld, W. J., Weerts, A. H., Visser, M., Buitink, J., Imhoff, R. O., Boisgontier,
H., Bouaziz, L., Eilander, D., Hegnauer, M., ten Velden, C., and Russell, B., 2024.
Wflow_sbm v0.7.3, a spatially distributed hydrological model: from global data to local
applications. Geosci. Model Dev., 17, 3199–3234. <https://doi.org/10.5194/gmd-17-3199-2024>.

To cite a specific software version of wflow please use the DOI provided in the Zenodo badge
above, that points to the latest release.