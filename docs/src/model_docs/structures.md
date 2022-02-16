# Julia structures

## Model

 Below the composite type that represents all different aspects of a `Wflow.Model`, such as
 the network, parameters, clock, model type, configuration and input and output.

```julia
struct Model{N,L,V,R,W,T}
    config::Config  # all configuration options
    network::N      # connectivity information, directed graph
    lateral::L      # lateral model that holds lateral state, moves along network
    vertical::V     # vertical model that holds vertical state, independent of each other
    clock::Clock    # to keep track of simulation time
    reader::R       # provides the model with dynamic input
    writer::W       # writes model output
    type::T         # model type
end
```

The `lateral` field of the `struct Model` can contain different lateral concepts. For each
Wflow model these different lateral concepts are mapped through the use of a `NamedTuple`.
The `vertical` field of the `struct Model` always contains one vertical concept, for example
the SBM or HBV vertical concept.

Below an example how lateral concepts are mapped for the SBM model through a `NamedTuple`:

```julia
(subsurface = ssf, land = olf, river = rf)
```

The `subsurface` part is mapped to the lateral subsurface flow kinematic wave concept, the
`land` part is mapped the overland flow kinematic wave concept and the `river` part is
mapped to the river flow kinematic wave concept. Knowledge of this specific mapping is
required to understand and correctly set input, output and state variables in the TOML
configuration file, see also [Config and TOML](@ref config_toml). This mapping is described in more
detail for each model in the section Models. Also the `struct` of each mapped concept is
provided, so one can check the internal variables in the code. These structs are defined as
a parametric composite type, with type parameters between curly braces after the `struct`
name. See also the next paragraph [Vertical and lateral models](@ref) for a more
detailed description.

## Vertical and lateral models
The different model concepts used in Wflow are defined as parametric [composite
types](https://docs.julialang.org/en/v1/manual/types/#Composite-Types). For example the
vertical `SBM` concept is defined as follows: `struct SBM{T,N,M}`. `T`, `N` and `M` between
curly braces after the `struct` name refer to type parameters, for more information about
type parameters you can check out [Type
parameters](https://docs.julialang.org/en/v1/manual/types/#man-parametric-composite-types).
Since these parameters can be of any type, it is possible to declare an unlimited number of
composite types. The type parameters are used to set the type of `struct` fields, below an
example with a part of the `SBM` struct:

```julia
@get_units @with_kw struct SBM{T,N,M}
    # Model time step [s]
    Δt::T | "s"
    # Maximum number of soil layers
    maxlayers::Int | "-"
    # number of cells
    n::Int | "-"
    # Number of soil layers
    nlayers::Vector{Int} | "-"
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int} | "-"
    # Fraction of river [-]
    riverfrac::Vector{T} | "-"
    # Saturated water content (porosity) [mm mm⁻¹]
    θₛ::Vector{T} | "mm mm-1"
    # Residual water content [mm mm⁻¹]
    θᵣ::Vector{T} | "mm mm-1"
    # Vertical hydraulic conductivity [mm Δt⁻¹] at soil surface
    kv₀::Vector{T} | "mm Δt-1"
    # Muliplication factor [-] applied to kv_z (vertical flow)
    kvfrac::Vector{SVector{N,T}} | "-"
```

The type parameter `T` is used in Wflow as a subtype of `AbstractFloat`, allowing to store
fields with a certain floating point precision (e.g. `Float64` or `Float32`) in a flexible
way. `N` refers to the maximum number of soil layers of the `SBM` soil column, and `M`
refers to the maximum number of soil layers + 1. See also part of the following instance of
`SBM`:

```julia
sbm = SBM{Float,maxlayers,maxlayers + 1}(
    Δt = tosecond(Δt),
    maxlayers = maxlayers,
    n = n,
```

For the other model concepts, we refer to the code to check these type parameters.