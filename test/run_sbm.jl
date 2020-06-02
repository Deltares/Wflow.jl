function update(model, toposort, n)
    @unpack lateral, vertical, network, clock = model

    Wflow.update_before_lateralflow(vertical)
    lateral.recharge[:] = vertical.recharge
    lateral.recharge .*= lateral.dl
    lateral.zi[:] = vertical.zi

    Wflow.update(lateral, network, toposort, n)

    Wflow.update_after_lateralflow(vertical, lateral.zi, lateral.exfiltwater)

    # update the clock
    clock.iteration += 1
    clock.time += clock.Δt

    return model
end

model = Wflow.initialize_sbm_model(staticmaps_moselle_path, leafarea_moselle_path)
toposort = Wflow.topological_sort_by_dfs(model.network)
n = length(toposort)
model = update(model, toposort, n)

@testset "first timestep" begin
    sbm = model.vertical
    @test model.clock.iteration == 2
    @test sbm.altitude[1] == 345.1470031738281
    @test sbm.θₛ[1] == 0.46367356181144714
    @test isnan(sbm.runoff[1])  # should probably be initialized to 0.0
    @test sbm.soilevap[1] == 0.08821308159314034
end

# run the second timestep
model = update(model, toposort, n)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.altitude[1] == 345.1470031738281
    @test sbm.θₛ[1] == 0.46367356181144714
    @test isnan(sbm.runoff[1])
    @test sbm.soilevap[1] == 0.17235604508244792
end

@testset "subsurface flow" begin
    ssf = model.lateral.ssf
    @test sum(ssf) ≈ 6.959580661699383e16
    @test ssf[toposort[1]] ≈ 4.392529226944353e11
    @test ssf[toposort[n-100]] ≈ 8.003673229321337e11
    @test ssf[sink] ≈ 6.92054650606041e11
end

using BenchmarkTools, Juno

benchmark = @benchmark update(model, toposort, n)
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
