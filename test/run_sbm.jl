tomlpath = joinpath(@__DIR__, "config.toml")
tomldir = dirname(tomlpath)
config = Wflow.Config(Wflow.parsefile(tomlpath))
output_path = joinpath(@__DIR__, "data", "output_run_sbm.nc")

model = Wflow.initialize_sbm_model(
    config,
    staticmaps_moselle_path,
    cyclic_moselle_path,
    forcing_moselle_path,
    instates_moselle_path,
    output_path,
)

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

    vcell = model.vertical[1]
    @test vcell isa NamedTuple
    @test isbits(vcell)
    @test vcell.tt ≈ 1.2999999523162842

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
    @test sbm.soilevap[1] == 0.00586565362774996
    @test sbm.snow[1] == 0.009696763863612956
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 7.005489495052358e16
    @test ssf[toposort_land[1]] ≈ 3.0449782003445332e13
    @test ssf[toposort_land[nl-100]] ≈ 7.87333555063647e11
    @test ssf[sink] ≈ 2.2998636307068414e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 6.1333024054146446
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[sink] ≈ 0.0
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 671.4362633660646
    @test q[4061] ≈ 0.011731014256037769
    @test q[5617] ≈ 1.0080691890749913
    @test q[toposort_river[end]] ≈ 0.0043612246315573805
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    inds = filter(i -> !ismissing(res[i]), 1:nr)
    @test res[inds[2]].outflow ≈ 0.2174998580279266
    @test res[inds[2]].inflow ≈ 0.004452979665488473
    @test res[inds[2]].volume ≈ 2.7745799188635208e7
    @test res[inds[2]].precipitation ≈ 3.0
    @test res[inds[2]].evaporation ≈ 4.0
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
