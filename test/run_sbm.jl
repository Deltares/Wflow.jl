function update(
    model,
    toposort_land,
    toposort_river,
    frac_toriver,
    index_river,
    nl,
    nr;
    masswasting = false,
)
    @unpack lateral, vertical, network, clock = model

    Wflow.update_forcing!(model)

    Wflow.update_until_snow(vertical)

    if masswasting
        snowflux_frac =
            min.(0.5, lateral.land.sl ./ 5.67) .* min.(1.0, vertical.snow ./ 10000.0)
        maxflux = snowflux_frac .* vertical.snow
        vertical.snow .=
            Wflow.accucapacityflux(network.land, toposort_land, vertical.snow, maxflux)
    end

    Wflow.update_until_recharge(vertical)
    lateral.subsurface.recharge .= vertical.recharge
    lateral.subsurface.recharge .*= lateral.subsurface.dl
    lateral.subsurface.zi .= vertical.zi

    Wflow.update(
        lateral.subsurface,
        network.land,
        toposort_land,
        frac_toriver,
        lateral.river.rivercells,
    )

    Wflow.update_after_lateralflow(
        vertical,
        lateral.subsurface.zi,
        lateral.subsurface.exfiltwater,
    )

    lateral.land.qlat .=
        (vertical.runoff .* vertical.xl .* vertical.yl .* 0.001) ./ 86400.0 ./
        lateral.land.dl

    Wflow.update(
        lateral.land,
        network.land,
        toposort_land,
        nl,
        frac_toriver = frac_toriver,
        river = lateral.river.rivercells,
        do_iter = true,
    )

    lateral.river.qlat .=
        (
            lateral.subsurface.to_river[index_river] ./ 1.0e9 ./ lateral.river.Δt .+
            lateral.land.to_river[index_river]
        ) ./ lateral.river.dl

    Wflow.update(lateral.river, network.river, toposort_river, nr, do_iter = true)

    Wflow.write_output(model, model.writer)

    # update the clock
    clock.iteration += 1
    clock.time += clock.Δt

    return model
end

tomlpath = joinpath(@__DIR__, "config.toml")
tomldir = dirname(tomlpath)
config = Wflow.Config(TOML.parsefile(tomlpath))
output_path = joinpath(@__DIR__, "data", "output_run_sbm.nc")

model = Wflow.initialize_sbm_model(
    config,
    staticmaps_moselle_path,
    cyclic_moselle_path,
    forcing_moselle_path,
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

model = update(
    model,
    toposort_land,
    toposort_river,
    frac_toriver,
    index_river,
    nl,
    nr,
    masswasting = true,
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
    @test sbm.soilevap[1] == 0.048453431759442424
    @test sbm.snow[1] == 0.0
end

# run the second timestep
model = update(
    model,
    toposort_land,
    toposort_river,
    frac_toriver,
    index_river,
    nl,
    nr,
    masswasting = true,
)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.altitude[1] == 643.5469970703125
    @test sbm.θₛ[1] == 0.48343977332115173
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.09563376762901096
    @test sbm.snow[1] == 0.0
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.8949450377727496e16
    @test ssf[toposort_land[1]] ≈ 3.0293145813424688e13
    @test ssf[toposort_land[nl-100]] ≈ 7.566873257040469e11
    @test ssf[sink] ≈ 2.1952330312819342e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 6.012832744996996
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[sink] ≈ 0.0
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 655.0918562695385
    @test q[4061] ≈ 0.011394640544977987
    @test q[5617] ≈ 0.986801851306535
    @test q[toposort_river[end]] ≈ 0.004621542892763028
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    inds = filter(i -> !isequal(res[i], nothing), 1:nr)
    @test res[inds[2]].outflow ≈ 0.2174998021586705
    @test res[inds[2]].inflow ≈ 0.07074779883466663
    @test res[inds[2]].volume ≈ 2.7385279244359117e7
    @test res[inds[2]].precipitation ≈ 3.0
    @test res[inds[2]].evaporation ≈ 4.0
end

using BenchmarkTools, Juno

benchmark = @benchmark update(
    model,
    toposort_land,
    toposort_river,
    frac_toriver,
    index_river,
    nl,
    nr,
    masswasting = true,
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
