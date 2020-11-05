# Model structure

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
