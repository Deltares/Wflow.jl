tomlpath = joinpath(@__DIR__, "config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)

toposort_land = Wflow.topological_sort_by_dfs(model.network.land)
toposort_river = Wflow.topological_sort_by_dfs(model.network.river)
nl = length(toposort_land)
nr = length(toposort_river)
index_river = filter(i -> !isequal(model.lateral.river.rivercells[i], 0), 1:nl)
frac_toriver = Wflow.fraction_runoff_toriver(
    model.network.land,
    index_river,
    model.lateral.subsurface.βₗ,
    nl,
)

model = Wflow.update(
    model,
    config,
    toposort_land,
    toposort_river,
    frac_toriver,
    index_river,
    nl,
    nr,
)

@testset "first timestep" begin
    sbm = model.vertical
    
    @test sbm.tt[1] ≈ 1.2999999523162842

    @test model.clock.iteration == 2

    @test sbm.altitude[1] == 643.5469970703125
    @test sbm.θₛ[1] == 0.48343977332115173
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.snow[1] == 0.6029994894902302
end

# run the second timestep
model = Wflow.update(
    model,
    config,
    toposort_land,
    toposort_river,
    frac_toriver,
    index_river,
    nl,
    nr,
)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.altitude[1] == 643.5469970703125
    @test sbm.θₛ[1] == 0.48343977332115173
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.005865653628210413
    @test sbm.snow[1] == 0.009696763863612956
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 7.005491155295535e16
    @test ssf[toposort_land[1]] ≈ 3.0449782003445332e13
    @test ssf[toposort_land[nl-100]] ≈ 7.87333555063647e11
    @test ssf[sink] ≈ 2.2998636307068414e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 6.132993832259066
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[sink] ≈ 0.0
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 671.4270822879581
    @test q[4061] ≈ 0.011731014256037769
    @test q[5617] ≈ 1.0080691483854958
    @test q[toposort_river[end]] ≈ 0.0043612246315573805
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    # the index of the second reservoir
    i = findall(res.is_res)[2]

    @test res.outflow[i] ≈ 0.2174998580279266
    @test res.inflow[i] ≈ 0.004452979665488473
    @test res.volume[i] ≈ 2.7745799188635208e7
    @test res.precipitation[i] ≈ 3.0
    @test res.evaporation[i] ≈ 4.0
end

benchmark = @benchmark Wflow.update(
    model,
    config,
    toposort_land,
    toposort_river,
    frac_toriver,
    index_river,
    nl,
    nr,
)

"Prints a benchmark results just like btime"
function print_benchmark(trialmin)
    trialtime = BenchmarkTools.time(trialmin)
    trialallocs = BenchmarkTools.allocs(trialmin)
    println(
        "  ",
        BenchmarkTools.prettytime(trialtime),
        " (",
        trialallocs,
        " allocation",
        trialallocs == 1 ? "" : "s",
        ": ",
        BenchmarkTools.prettymemory(BenchmarkTools.memory(trialmin)),
        ")",
    )
end

trialmin = BenchmarkTools.minimum(benchmark)

println("Model update")
print_benchmark(trialmin)
# @profview Wflow.update(model, config, toposort_land, toposort_river, frac_toriver, index_river, nl, nr)

Wflow.close_files(model)
