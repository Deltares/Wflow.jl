# ensure test data is present
datadir = joinpath(@__DIR__, "data")
isdir(datadir) || mkdir(datadir)
staticmaps_moselle_path = joinpath(datadir, "staticmaps-moselle.nc")
isfile(staticmaps_moselle_path) || download("https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/staticmaps.nc", staticmaps_moselle_path)
leafarea_moselle_path = joinpath(datadir, "lai_clim-moselle.nc")
isfile(leafarea_moselle_path) || download("https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/lai_clim.nc", leafarea_moselle_path)


function update(model, toposort, n)
    # update the vertical
    for i in eachindex(model.vertical)
        model.vertical[i] = Wflow.update(model.vertical[i])
        model.lateral.recharge[i] = model.vertical[i].recharge  * model.lateral.dl[i]
        model.lateral.zi[i] = model.vertical[i].zi
    end
    model.lateral = Wflow.update(model.lateral, model.network, toposort, n)
    # update the clock
    model.clock.iteration += 1
    model.clock.time += model.clock.Δt

    return model
end

model = Wflow.initialize_sbm_model(staticmaps_moselle_path, leafarea_moselle_path)
toposort = Wflow.topological_sort_by_dfs(model.network)
n = length(toposort)
model = update(model, toposort, n)
