# Model structure

 Below the composite type that represents all different aspects of a `Wflow.Model`, such as
 the network, parameters, clock, configuration and input and output.

```julia
struct Model{N,L,V,R,W}
    config::Config  # all configuration options
    network::N      # connectivity information, directed graph
    lateral::L      # lateral model that holds lateral state, moves along network
    vertical::V     # vertical model that holds vertical state, independent of each other
    clock::Clock    # to keep track of simulation time
    reader::R       # provides the model with dynamic input
    writer::W       # writes model output
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
configuration file, see also [Config and TOML](@ref). This mapping is described in more
detail for each model in the section Models. Also the `struct` of each mapped concept is
provided, so one can check the internal variables in the code.