using LightGraphs
using Dates

mutable struct Clock
    time::DateTime
    iteration::Int
    Δt::Second
end

struct Model{N,L,V,R,W}
    network::N  # connectivity information, directed graph
    lateral::L  # lateral model that holds lateral state, moves along network
    vertical::Vector{V}  # vertical model that holds vertical state, independent of each other
    clock::Clock  # to keep track of simulation time
    reader::R  # provides the model with dynamic input
    writer::W  # writes model output
end
