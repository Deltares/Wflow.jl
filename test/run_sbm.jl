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

    # increases time from 115ms to 3.6s and allocations from 0 to 150MB
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
    lateral.subsurface.recharge[:] = vertical.recharge
    lateral.subsurface.recharge .*= lateral.subsurface.dl
    lateral.subsurface.zi[:] = vertical.zi

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

    lateral.land.qlat[:] =
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

    lateral.river.qlat[:] =
        (
            lateral.subsurface.to_river[index_river] ./ 1.0e9 ./ lateral.river.Δt .+
            lateral.land.to_river[index_river]
        ) ./ lateral.river.dl

    Wflow.update(lateral.river, network.river, toposort_river, nr, do_iter = true)

    # update the clock
    clock.iteration += 1
    clock.time += clock.Δt

    return model
end

writer = nothing  # TODO use a CSV writer, RowWriter?

model = Wflow.initialize_sbm_model(
    staticmaps_moselle_path,
    leafarea_moselle_path,
    forcing_moselle_path,
    writer,
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
    @test sbm.soilevap[1] == 0.07996413903665732
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
    @test sbm.soilevap[1] == 0.15658817917471723
    @test sbm.snow[1] == 0.0
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.910372076543659e16
    @test ssf[toposort_land[1]] ≈ 3.0271045626104652e13
    @test ssf[toposort_land[nl-100]] ≈ 7.555413711433021e11
    @test ssf[sink] ≈ 2.177681570044858e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 6.044737022735466
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[sink] ≈ 0.0
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 656.6441748952489
    @test q[4061] ≈ 0.011400418905152046
    @test q[5617] ≈ 0.9858388071125663
    @test q[toposort_river[end]] ≈ 0.004344447476624051
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    inds = filter(i -> !isequal(res[i], nothing), 1:nr)
    @test res[inds[2]].outflow ≈ 0.21749979713718748
    @test res[inds[2]].inflow ≈ 0.051269540649685884
    @test res[inds[2]].volume ≈ 2.7357871245776203e7
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
# @time update(model, toposort, n)
# @btime update(model, toposort, n)
# @profiler foreach(x -> update(model, toposort, n), 1:10)  # run a few times for more accuracy

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
trialallocs = BenchmarkTools.allocs(trialmin)
trialtime = BenchmarkTools.time(trialmin)

if VERSION >= v"1.5"
    # views don't allocate anymore in julia 1.5
    @test trialallocs == 0
end
println("Model update")
print_benchmark(trialmin)
