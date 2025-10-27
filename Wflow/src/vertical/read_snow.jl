function NewSnowModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    n = length(indices)
    cache = NewSnowCache(; n)
    properties = SnowHbvParameters(dataset, config, indices, dt)
    p = NewSnowParameters(; cache, properties)
    state_prototype = ComponentVector(;
        snow_storage = 0.0,
        snow_water = 0.0,
        cumulative_snow_melt = 0.0,
        cumulative_runoff = 0.0,
    )
    integrators = [WflowIntegrator(copy(state_prototype)) for _ in 1:n]
    return NewSnowModel(; p, integrators)
end